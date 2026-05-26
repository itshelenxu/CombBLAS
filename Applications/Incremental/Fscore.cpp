#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <string.h>
#include <assert.h>
#include <omp.h>
#include <map>
using namespace std;

int significantClusterSize = 0;

template <typename T>
vector<size_t> sortIndices(const vector<T> &v) {
    
    // initialize original index locations
    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);
    
    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
         [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});
    
    return idx;
}



double Fscore(string hipmclFile, string incmclFile, int base)
{
    ifstream hipmclStream (hipmclFile);
    ifstream incmclStream (incmclFile);
    string line;
    int64_t item, clustID = 0;
    int64_t numIncMCLClusters, numHipMCLClusters, numSigHipMCLClusters = 0;
    int64_t nproteins = 0;
    int64_t nproteinsTarget = 0;
    // item 0-based
    std::map<int64_t, int64_t> vtxMap;
    std::map<int64_t, int64_t> vtxMapR;
    int64_t vtxCount=0;
    
    vector<vector<int64_t>> hipmclClusters;
    
    if (hipmclStream.is_open())
    {
        while ( getline (hipmclStream,line) )
        {
            istringstream iss(line);
            if(clustID >= hipmclClusters.size())
            {
                hipmclClusters.resize(clustID * 2 + 1);
            }
            while ( iss >> item)
            {
                //cout << item << endl;
                if(base==1) item--;
                vtxMap[item] = vtxCount;
                vtxMapR[vtxCount] = item;
                hipmclClusters[clustID].push_back(vtxCount);
                vtxCount = vtxCount + 1;
            }
            
            nproteins += hipmclClusters[clustID].size();
            if (hipmclClusters[clustID].size() >= significantClusterSize){
                nproteinsTarget += hipmclClusters[clustID].size(); // Target only those proteins who belong to significant clusters in HipMCL
                numSigHipMCLClusters++;
            }
            clustID++;
        }
        hipmclStream.close();
        numHipMCLClusters = clustID;
    }
    else
    {
        cout << "Unable to open " << hipmclFile << endl;
        return -1;
    }
    std::vector<int64_t> hipmclClusterAsn(nproteins, -1);
    for(int64_t i = 0; i < hipmclClusters.size(); i++){
        for(int64_t j =0 ; j < hipmclClusters[i].size(); j++) hipmclClusterAsn[hipmclClusters[i][j]] = i;
    }
    
    cout << "Number of clusters from HipMCL: " << numHipMCLClusters << endl;
    cout << "Number of significant clusters from HipMCL: " << numSigHipMCLClusters << endl;
    cout << "Number of target proteins: " << nproteinsTarget << endl;
    vector<int64_t> incmclClusterAsn(nproteins, -1);
    clustID = 0;
    if (incmclStream.is_open())
    {
        while ( getline (incmclStream,line) )
        {
            istringstream iss(line);
            while ( iss >> item)
            {
                if(base==1) item--;
                //int64_t hipmclclusID = hipmclClustAsn[];
                int64_t itemMapped = vtxMap.find(item)->second;
                if(hipmclClusters[hipmclClusterAsn[itemMapped]].size() >= significantClusterSize) incmclClusterAsn[itemMapped] = clustID;
                //if(vtxMap.find(item)->second >= nproteins)
                //{   
                    //cout << item << "-" << nproteins << endl;
                    //cout << "The number of vertices in HipMCL and IncMCL outputs does not match. \nExiting.............." << endl;
                    //exit(1);
                //}
            }
            clustID++;
        }
        incmclStream.close();
        numIncMCLClusters = clustID;
    }
    else
    {
        cout << "Unable to open " << incmclFile << endl;
        return -1;
    }
    
    // this will account for isolated vertices. These vertices are assigned to their own clusters 
    for(int i=0; i<nproteins; i++)
    {
        if(incmclClusterAsn[i]==-1) incmclClusterAsn[i]=numIncMCLClusters++;
    }
    //cout << "Number of clusters from IncMCL: " << numIncMCLClusters << endl;
    
    vector<int64_t> clusterSizes2(numIncMCLClusters,0);
    for(int i=0; i<nproteins; i++)
    {
        clusterSizes2[incmclClusterAsn[i]]++;
    }
    
    
    vector<double> F(numHipMCLClusters);
    vector<int64_t> nisect(numIncMCLClusters); // number of items in the intersection
    vector<int64_t> isect(numIncMCLClusters); // clusters with nonzero intersection with the current cluster
    int mismatch = 0;

    for(int i=0; i<numHipMCLClusters; i++)
    {
        if(hipmclClusters[i].size() >= significantClusterSize){
            int64_t isectCount = 0;
            fill(nisect.begin(), nisect.end(), 0);
            for(int j=0; j<hipmclClusters[i].size(); j++)
            {
                int64_t item1 = hipmclClusters[i][j];
                int64_t c2 = incmclClusterAsn[item1];
                if(nisect[c2]==0) isect[isectCount++] = c2;
                nisect[c2] ++;
            }
            auto maxoverlap = max_element(nisect.begin(), nisect.end());
            int64_t c2_max = distance(nisect.begin(), maxoverlap);
            if(*maxoverlap!=hipmclClusters[i].size() || *maxoverlap!=clusterSizes2[c2_max])
            {
                cout << "Mismatch# " << mismatch++ << ":: HipMCL Cluster: "<< i << " Size: " << hipmclClusters[i].size() << " IncMCL Cluster: "<< c2_max <<" Size: " << clusterSizes2[c2_max] << endl;
                
            }
            double Fi = 0;
            for(int j=0; j<isectCount; j++)
            {
                int64_t c2 = isect[j];
                double precision = (double) nisect[c2] / clusterSizes2[c2];
                double recall = (double) nisect[c2] / hipmclClusters[i].size();
                double Fij = 2 * precision * recall / (precision + recall);
                Fi = max(Fi, Fij);
            }
            
            Fi = Fi * hipmclClusters[i].size()/ nproteinsTarget;
            F[i] = Fi;
        }
    }
    return accumulate(F.begin(), F.end(), 0.0);
}


int main(int argc, char* argv[])
{
    
    string ifilename1 = "";
    string ifilename2 = "";
    int base = 0;
    int threshold = 0;
  
    
    if(argc != 9)
    {
        cout << "Usage: ./fscore -M1 <MCLOut> -M2 <HipMCLOut> -base <Base of vertices (same as HipMCL) 1 or 0> -threshold <sizeof HipMCL cluster>\n";
        cout << "Example: ./fscore -M1 input1.txt -M2 input2.txt -base 0 -threshold 0" << endl;
        return -1;
    }
    
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i],"-M1")==0)
        {
            ifilename1 = string(argv[i+1]);
            printf("\nMCL output file: %s",ifilename1.c_str());
        }
        else if (strcmp(argv[i],"-M2")==0)
        {
            ifilename2 = string(argv[i+1]);
            printf("\nHipMCL output file: %s",ifilename2.c_str());
        }
        else if (strcmp(argv[i],"-base")==0)
        {
            base = atoi(argv[i+1]);
            printf("\nbase: %d",base);
        }
        else if (strcmp(argv[i],"-threshold")==0)
        {
            threshold = atoi(argv[i+1]);
            printf("\nthreshold of significant cluster size: %d",threshold);
            significantClusterSize = threshold;
        }
    }
    printf("\n");
    double F = Fscore(ifilename1, ifilename2, base);
    cout << "F score: " << F << endl;
    return 0;
}
