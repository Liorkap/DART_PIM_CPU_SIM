
#include "manager.h"

/*============================================= Heap ===============================================*/
/*==================================================================================================*/

void Heap::heapify(size_t index) {
    size_t size = elements.size();
    size_t largest = index;
    size_t left = 2 * index + 1;
    size_t right = 2 * index + 2;

    if (left < size && (elements[left]->kmerSeq < elements[largest]->kmerSeq)) {
        largest = left;
    }

    if (right < size && (elements[right]->kmerSeq < elements[largest]->kmerSeq)) {
        largest = right;
    }

    if (largest != index) {
        std::swap(elements[index], elements[largest]);
        heapify(largest);
    }
}

Heap::Heap(const std::vector<Kmer*>& values) : elements(values) {
    std::make_heap(elements.begin(), elements.end(), compareMinimizers);
}

void Heap::insert(Kmer* value) {
    elements.push_back(value);
    size_t index = elements.size() - 1;

    while (index > 0) {
        size_t parent = (index - 1) / 2;
        if (elements[index]->kmerSeq > elements[parent]->kmerSeq) {
            break;
        }

        std::swap(elements[index], elements[parent]);
        index = parent;
    }
    heapify(0);
}

Kmer* Heap::minElement() {
    return elements[0];
}

bool Heap::compareMinimizers(const Kmer* a, const Kmer* b) {
    return b->kmerSeq < a->kmerSeq;
}

/*==================================================================================================*/
/*==================================================================================================*/









/*============================================= Kmer ===============================================*/
/*==================================================================================================*/

Kmer::Kmer(int pos, const string &kmerSeq) : position(pos), kmerSeq(kmerSeq) {}
Kmer::Kmer(const Kmer &kmer): position(kmer.position), kmerSeq(kmer.kmerSeq) {}

/*==================================================================================================*/
/*==================================================================================================*/








/*======================================= RefGenomeMinimizer =======================================*/
/*==================================================================================================*/

RefGenomeMinimizer::RefGenomeMinimizer(Kmer minimizer, string refSegment): minimizer(minimizer), refSegment(refSegment) {
    refSegmentPosition = minimizer.position - (READ_LENGTH- KMER_LENGTH + ERROR_THRESHOLD);
}


int RefGenomeMinimizer::getWFSeq(int readMinimizerPosition, string* refSubSeq){
    int refSubSeqPosition = minimizer.position - readMinimizerPosition; // The location of the ref sub sequence relative to the whole genome;
    refSubSeqPosition -= refSegmentPosition; // The location of the ref sub sequence relative to the refSegment;
    refSubSeqPosition -= ERROR_THRESHOLD; // Adding the error threshhold
    *refSubSeq = refSegment.substr(refSubSeqPosition, REF_SUB_SEQUENCE_LENGTH);
    return refSubSeqPosition + refSegmentPosition; // location in refGenome
}

void RefGenomeMinimizer::print(){
    std::cout << "minimizer: " << minimizer.kmerSeq << endl;
    std::cout << "minimizer position (in reference genome): " << minimizer.position << endl;
    std::cout << "reference segment: " << refSegment << endl;
}

/*==================================================================================================*/
/*==================================================================================================*/






/*========================================= ReadMinimizer ==========================================*/
/*==================================================================================================*/

ReadMinimizer::ReadMinimizer(Kmer minimizer) : minimizer(minimizer) {}

void ReadMinimizer::print(){
    std::cout << "---------------" << endl;
    std::cout << "minimizer: " << minimizer.kmerSeq << endl;
    std::cout << "minimizer position (in read): " << minimizer.position << endl;
    std::cout << "read potential location: " << readPotentialLocation << endl;
    std::cout << "score: " << score << endl;
    std::cout << "---------------" << endl;
}

/*==================================================================================================*/
/*==================================================================================================*/







/*============================================= Read ===============================================*/
/*==================================================================================================*/

Read::Read(string readSeq) {
    vector<Kmer*> readMinimizers = createMinimizers(readSeq);
    this->seq = readSeq;

    for(Kmer* minimizer : readMinimizers){
        ReadMinimizer readMinimizer(*minimizer);
        this->minimizers.push_back(readMinimizer);
    }
}

std::vector<Kmer*> Read::createMinimizers(const string &seq) {
    int i;
    std::vector<Kmer*> heapKmers;
    std::vector<Kmer*> kmersInOrder;
    Kmer* minMinimizer;
    std::vector<Kmer*> outMinimizers;

    for(i = 0; (i < seq.size()) && (i < WINDOW_SIZE); i++) {
        std::string subString = seq.substr(i, KMER_LENGTH);

        minMinimizer = new Kmer(i, subString);
        heapKmers.push_back(minMinimizer);
        kmersInOrder.push_back(minMinimizer);
    }
    Heap heapMin(heapKmers);

    for(; (i <= seq.size() - KMER_LENGTH + 1); i++) {
        auto currentMin = heapMin.minElement();
        if (outMinimizers.empty()) {
            outMinimizers.push_back(new Kmer(currentMin->position, currentMin->kmerSeq));
        }
        else if(currentMin->position != outMinimizers.back()->position) {
            outMinimizers.push_back(new Kmer(currentMin->position, currentMin->kmerSeq));
        }
        std::string subString = seq.substr(i, KMER_LENGTH);
        auto newMinimizer = new Kmer(i, subString);
        kmersInOrder.push_back(newMinimizer);
        auto elem1 = kmersInOrder[kmersInOrder.size() - WINDOW_SIZE - 1];
        auto iteratorToDelete = std::find(heapMin.elements.begin(), heapMin.elements.end(), elem1);
        heapMin.elements.erase(iteratorToDelete);

        heapMin.insert(newMinimizer);
    }
    return outMinimizers;
}

void Read::print(){
    std::cout << "--------------Printing Read Data--------------" << endl;
    std::cout << "Read Sequence: " << seq << endl;
    std::cout << "Minimizers:" << endl;
    for(ReadMinimizer minimizer : minimizers){
        minimizer.print();
    }
    std::cout << "----------------------------------------------" << endl;
}

/*==================================================================================================*/
/*==================================================================================================*/








/*========================================== PendingJob ============================================*/
/*==================================================================================================*/

PendingJob::PendingJob(string readSeq, string refSeq, int* score) : readSeq(readSeq), refSeq(refSeq), score(score) {}

/*==================================================================================================*/
/*==================================================================================================*/






/*========================================== PimPacket =============================================*/
/*==================================================================================================*/

PimPacket::PimPacket(string readSeq) : readSeq(readSeq) {}

/*==================================================================================================*/
/*==================================================================================================*/









/*============================================= Manager ============================================*/
/*==================================================================================================*/
Manager::Manager(CPUMinimizers CPUMins, vector<Read> reads){
    this->CPUMins = CPUMins;
    this->reads = reads;
    this->numRunningJobs = 0;
}

void Manager::handleReads(){
    // go over all reads
    for(Read& read : reads) {
        PimPacket pimData(read.seq);
        // go over all minimizers of a read
        for(ReadMinimizer& readMinimizer : read.minimizers){
            bool foundMinmizer = false;
            // go over the CPU minimizers and check if this minimizer is there
            for(int i = 0; i < CPUMins.size(); i++){
                auto refMinimizer = CPUMins[i];
                if(refMinimizer.minimizer.kmerSeq == readMinimizer.minimizer.kmerSeq){ // checking if its a CPU minimizer
                    string refSeq;
                    // get the sub reference seqment to send to WF
                    readMinimizer.readPotentialLocation = refMinimizer.getWFSeq(readMinimizer.minimizer.position, &refSeq);
                    // if there is a slot for a new WF job, send it
                    if(numRunningJobs < MAX_RUNNING_JOBS){

                        runningJobsMtx.lock();
                        numRunningJobs++;
                        runningJobsMtx.unlock();

                        thread WFJob(&Manager::wagnerFischerAffineGap, this, read.seq, refSeq, &readMinimizer.score, true, 1, 1, 1);
                        WFJob.detach();
                    }
                    // if not, add to pending jobs and call the pending jobs handler
                    else{
                        std::cout << "reached max running jobs" << endl;
                        PendingJob currRead(read.seq, refSeq, &readMinimizer.score);
                        pendingJobsForWF.push_back(currRead);
                        if(pendingJobsForWF.size() == 1){
                            handlePendingReads();
                        }
                    }
                    foundMinmizer = true;
                }
            }
            // if this minimizer is not a CPU minimizer add it to tthe pim packet
            if(!foundMinmizer){
                readMinimizer.score = -1;
                pimData.minimizers.push_back(readMinimizer.minimizer);
            }
        }
        // if some minimizers for DART-PIM are found send it to DART-PIM
        if(pimData.minimizers.size() > 0){
            //TODO: send the read and pim minimizers to pim if possible
        }
    }
}

void Manager::handlePendingReads(){
    // while we have pending WF jobs
    while(pendingJobsForWF.size() > 0){
        //check if there is a free slot to run a job
        if(numRunningJobs < MAX_RUNNING_JOBS){
            // get the next job
            PendingJob nextRead = pendingJobsForWF.front();
            pendingJobsForWF.pop_front();

            runningJobsMtx.lock();
            numRunningJobs++;
            runningJobsMtx.unlock();

            thread WFJob(&Manager::wagnerFischerAffineGap, this, nextRead.readSeq,  nextRead.refSeq, nextRead.score,
                         true, 1, 1, 1);
            WFJob.detach();
        }
    }
}

int Manager::wagnerFischerAffineGap(const string& S1, const string& S2, int* score,  bool backtraching, int wop, int wex, int wsub) {
    int n = S1.size();
    int m = S2.size();
    vector<int> contenders;
    int minCon = 0;
    string seq1_align = "";
    string seq2_align = "";

    // Initialize matrices D, M1, and M2 with appropriate dimensions
    vector<vector<int>> D(n + 1, vector<int>(m + 1, 0));
    vector<vector<int>> M1(n + 1, vector<int>(m + 1, 0));
    vector<vector<int>> M2(n + 1, vector<int>(m + 1, 0));

    // Initialize matrices with appropriate values
    for (int i = 1; i <= n; ++i) {
        D[i][0] = i * wex;
        M1[i][0] = i * wex;
    }
    for (int j = 1; j <= m; ++j) {
        D[0][j] = j * wex;
        M2[0][j] = j * wex;
    }

    // Fill in the matrices using dynamic programming
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            M1[i][j] = min(M1[i - 1][j] + wex, D[i - 1][j] + wop + wex);
            M2[i][j] = min(M2[i][j - 1] + wex, D[i][j - 1] + wop + wex);

            if (S1[i - 1] == S2[j - 1]) {
                D[i][j] = D[i - 1][j - 1];
            } else {
                D[i][j] = min({M1[i][j], M2[i][j], D[i - 1][j - 1] + wsub});
            }
        }
    }
    // Return the optimal alignment score
    *score = D[n][m];

    // find the alignment of the read according to sub reference sequence
    if(backtraching) {
        int i = n, j = m;
        while (i > 0 && j > 0 && D[i][j] > 0) {
            contenders = {D[i - 1][j - 1], D[i - 1][j], D[i][j - 1]};
            minCon = *(std::min_element(contenders.begin(), contenders.end()));

            if (minCon == contenders[0]) {
                if (S1[i - 1] != S2[j - 1]) {
                    seq1_align += S2[i - 1];
                    seq2_align += S2[j - 1];
                } else {
                    seq1_align += S1[i - 1];
                    seq2_align += S2[j - 1];
                }
                i -= 1;
                j -= 1;
            } else if (minCon == contenders[1]) {
                seq1_align += S1[i - 1];
                seq2_align += "-";
                i -= 1;
            } else {
                seq1_align += "-";
                seq2_align += S2[j - 1];
                j -= 1;
            }

        }
        std::reverse(seq1_align.begin(), seq1_align.end());
        std::reverse(seq2_align.begin(), seq2_align.end());
    }

    runningJobsMtx.lock();
    numRunningJobs--;
    runningJobsMtx.unlock();

    return *score;
}

void Manager::printCPUMinimizers(){
    int i = 0;
    std::cout << "--------------Printing CPU Minimizers--------------" << endl;
    for(RefGenomeMinimizer CPUMin : CPUMins){
        std::cout << "CPU Minimizer " << i << ":" << endl;
        CPUMin.print();
        i++;
    }
    std::cout << "---------------------------------------------------" << endl;
}

void Manager::printReads(){
    int i = 0;
    std::cout << "--------------Printing Reads--------------" << endl;
    for(Read read : reads){
        std::cout << "Read " << i << ":" << endl;
        read.print();
        i++;
    }
    std::cout << "------------------------------------------" << endl;
}

/*==================================================================================================*/
/*==================================================================================================*/








/*============================================= Main ===============================================*/
/*==================================================================================================*/

string randomizeSeq(int len){

    // String to store the random characters
    string result;

    // Generate random characters and append them to the string
    for (int i = 0; i < len; ++i) {
        char randomChar = '0' + (rand() % 4); // Random character between '0' and '3'
        result.push_back(randomChar);
    }

    return result;
}

vector<Read> getRandomReads(int numOfReads){
    vector<Read> reads;
    for(int i = 0 ; i < numOfReads; i++){
        Read read(randomizeSeq(READ_LENGTH));
        reads.push_back(read);
    }
    return reads;
}

string buildRandomRefSegment(string minimizer){
    string start = randomizeSeq(READ_LENGTH - KMER_LENGTH + ERROR_THRESHOLD);
    string end   = randomizeSeq(READ_LENGTH - KMER_LENGTH + ERROR_THRESHOLD);
    return start + minimizer + end;
}

CPUMinimizers getRandomCPUMinimizers(vector<Read> reads){
    CPUMinimizers randCPUMinimizers;


    // Generate a random number of reads to select
    int numReads = rand() % reads.size() +1;
    //std::cout << "num reads = " << numReads << endl;

    // Shuffle the indices to select minimizers randomly
    vector<int> readsIndices(reads.size());
    for(int i = 0; i < readsIndices.size(); i++){
        readsIndices[i] = i;
    }
    random_shuffle(readsIndices.begin(), readsIndices.end());

    for(int i = 0; i < numReads; i++){
        vector<ReadMinimizer> minimizers = reads[readsIndices[i]].minimizers;
        //std::cout << "chosen read = " << reads[readsIndices[i]].seq << endl;
        // Generate a random number of minimizers to select
        int numMinimizers = rand() % minimizers.size() +1; // Random number between 0 and minimizers.size()
        //std::cout << "num mins = " << numMinimizers << endl;

        // Shuffle the indices to select minimizers randomly
        vector<int> indices(minimizers.size());
        for (size_t i = 0; i < indices.size(); i++) {
            indices[i] = i;
        }
        random_shuffle(indices.begin(), indices.end());

        // Select random elements from the original vector based on the shuffled indices
        for (int i = 0; i < numMinimizers; i++) {
            string refSeg = buildRandomRefSegment(minimizers[indices[i]].minimizer.kmerSeq);
            Kmer readMinimizer = minimizers[indices[i]].minimizer;
            int refGenomeMinimizerMinLocation = READ_LENGTH - KMER_LENGTH + ERROR_THRESHOLD;
            int refGenomeMinimizerMaxLocation = REF_GENOME_LENGTH - (READ_LENGTH + ERROR_THRESHOLD);
            readMinimizer.position = refGenomeMinimizerMinLocation + rand() % (refGenomeMinimizerMaxLocation - refGenomeMinimizerMinLocation + 1);
            RefGenomeMinimizer refMinimizer(readMinimizer, refSeg);
            randCPUMinimizers.push_back(refMinimizer);
        }
    }

    return randCPUMinimizers;
}


void print_help(){
    std::cout << "------------------------------------------------ HELP ------------------------------------------------------" << endl;
    std::cout << "Two running modes:" << endl;
    std::cout << "   ./DART_PIM -rand                                                 :   will randomize reads and CPU minimizers" << endl;
    std::cout << "   ./DART_PIM -reads /path/to/reads/file -mins /path/to/mins/file   :   will run on cpecified reads and mins" << endl;
    std::cout << "------------------------------------------------------------------------------------------------------------" << endl;
}

void convertSeq2Nums(string& seq){
    for(int i = 0; i < seq.size(); i++){
        switch(seq[i]){
            case 'A':
                seq[i] = '0' + A;
                break;
            case 'C':
                seq[i] = '0' + C;
                break;
            case 'G':
                seq[i] = '0' + G;
                break;
            case 'T':
                seq[i] = '0' + T;
                break;
            default:
                std:cout << "Can't convert seq element " << seq[i] << "to number representasion" << endl;
        }
    }
}

void getReadsFromFile(ifstream& readsFile, vector<Read>& reads){
    string line;
    //skip first line
    if(!getline(readsFile, line)){
        std::cout << "MSG: Reads file is empty." << line << endl;
        return;
    }
    while(getline(readsFile, line)){
        convertSeq2Nums(line);
        Read read(line);
        reads.push_back(read);



        //skip two lines
        for(int i = 0; i < 3; i++){
            if(!getline(readsFile, line)){
                break;
            }
        }
    }
}

void getCPUMinsFromFile(ifstream& minsFile, CPUMinimizers& CPUMins){
    //TODO
}


int main(int argc, char* argv[]) {
    vector<Read> reads;
    CPUMinimizers CPUMins;
    ifstream readsFile;
    ifstream minsFile;
    bool readsFileOpen = false;
    bool minsFileOpen = false;
    int numOfReads = 100; //relevant to the rand running option
    argc = 2;
    (argv[1])= "-rand";

    if(argc == 2 && string(argv[1]) == "-rand"){


        srand(time(0));

        reads = getRandomReads(numOfReads);
        CPUMins = getRandomCPUMinimizers(reads);

    }
    else if(argc == 5 && string(argv[1]) == "-reads" && string(argv[3]) == "-mins"){
        readsFile = ifstream(argv[2]);
        readsFileOpen = readsFile.is_open();
        if(!readsFileOpen){
            std::cout << "ERROR: Can't open file " << string(argv[2]) << endl;
            return 1;
        }

        minsFile = ifstream(argv[4]);
        minsFileOpen = minsFile.is_open();
        if(!minsFileOpen){
            std::cout << "ERROR: Can't open file " << string(argv[4]) << endl;
            return 1;
        }

        getReadsFromFile(readsFile, reads);

        getCPUMinsFromFile(minsFile, CPUMins);
    }
    else if(argc == 2 && string(argv[1]) == "-help"){
        print_help();
        return 0;
    }
    else{
        std::cout << "ERROR: Invalid Arguments. See help: (run with -help)" << endl;
        print_help();
        return 1;
    }


    Manager manager(CPUMins, reads);

    manager.printCPUMinimizers();

    manager.handleReads();

    manager.printReads();



    if(readsFileOpen){
        readsFile.close();
    }

    if(minsFileOpen){
        minsFile.close();
    }

    return 0;
}

/*==================================================================================================*/
/*==================================================================================================*/