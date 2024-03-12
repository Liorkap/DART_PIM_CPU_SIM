
#define MAX_RUNNING_JOBS 3 // Max threads CPU can handle
// Nucleotides
#define A 0
#define C 1
#define G 2
#define T 3

#define KMER_LENGTH                12 
#define WINDOW_SIZE                4 //TODO: What window size does Rotem use?
#define READ_LENGTH                150
#define ERROR_THRESHOLD            3
#define REF_SUB_SEQUENCE_LENGTH    READ_LENGTH + 2 * ERROR_THRESHOLD
#define REF_GENOME_LENGTH          3000000000
#define UINT_MAX (INT_MAX * 2U + 1U)

#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <deque>
#include <map>
#include <thread>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <mutex>
#include <sstream>


using namespace std;

class Kmer {
public:
    int position;
    string kmerSeq;

    Kmer(int, const string &);
    Kmer(const Kmer&);
};

class Heap {
public:
    std::vector<Kmer*> elements;

    void heapify(size_t);
    Heap(const std::vector<Kmer*>&);
    void insert(Kmer*);
    void eraseLastMin();
    Kmer* minElement();
    static bool compareMinimizers(const Kmer*, const Kmer*);
};

//This is a row in the CPU minimizers list
class RefGenomeMinimizer {
public:
    Kmer   minimizer;  // Position is relative to the whole reference genome.     
    string refSegment; // Reference segment, the minimizer is in the middle (sort of).
    int    refSegmentPosition; // The position of the refSegment relative to the whole genome.

    RefGenomeMinimizer(Kmer, string);
    int getWFSeq(int, string*); // Extracts the relevant subRefSegment to send to WF. returns the potential location in refGenome. 
    void print();
};

// This is the CPU minimizers list type
typedef vector<RefGenomeMinimizer> CPUMinimizers;

// This is a type for a minimizer of a read
// Eventually it will get a score and a potenatial_location in the ref genome(either from CPU or DART-PIM)
class ReadMinimizer{
public:
    Kmer minimizer; // Position is relative to the read. 
    int score;
    int readPotentialLocation;

    ReadMinimizer(Kmer);
    void print();
};

// This is a type for a read. it hase a minimizers list and a location (will be updated eventualy based on the scores of the minimizers)
class Read{
public:
    string                seq;
    vector<ReadMinimizer> minimizers;
    int                   location;  // This is the actual result (location in the whole refrence genome).

    Read(string);
    vector<Kmer*> createMinimizers(const string &seq);
    std::vector<Kmer> findMinimizers(string s);
    void print();
};

// This is a type for a WF job that couldn't be scheduled because max number of threads were running.
// We store here the relevant data and schedule the job later when possible.
class PendingJob{
public:
    string readSeq;
    string refSeq;
    int*   score;

    PendingJob(string, string, int*);
};

// This is a type for the packet sent to the DART-PIM for calculation.
class PimPacket{
public:
    string       readSeq;
    vector<Kmer> minimizers;

    PimPacket(string);
};


//This is the manager type. 
class Manager{
public:

    CPUMinimizers       CPUMins; // List of spu minimizers
    vector<Read>        reads; // Reads to handle
    deque<PendingJob>   pendingJobsForWF; // These are the pending WF jobs that couldn't be scheduled and will be scheduled when possible. (too many parallel threads)
    int                 numRunningJobs; // Holds the current number of WF running threads
    mutex               runningJobsMtx; // The WF function reduces the number of running threads when its done. This is a shared variable across all WF jobs, protecting it with mutex. 

    Manager(CPUMinimizers, vector<Read>);
    void handleReads(); // handling the reads
    void handlePendingReads(); // handling the pending WF jobs
    int wagnerFischerAffineGap(const string& S1, const string& S2, int* score, bool backtracking, int wop=1, int wex=1, int wsub=1);
    void printReads();
    void printCPUMinimizers();
};

