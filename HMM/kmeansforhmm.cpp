#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <conio.h>
#include <sstream>
#include "config.h" //include configuration file that helps in global settings
#include <list>
#include <iterator>
#include <sstream>
#include <iomanip> //to print double with desired precision
#include <time.h> //use timer function to generate actual random indices

using namespace std;

long cepcoecount=0; //count the number of cepstral coefficients
double distortion=0,olddistortion=0; //distortion-->current distortion
int iterations=0;
string line;
double w[cepsize];//to store tokhura weights given by sir
int clustrep[cbsize]={0}; //word represented by the cluster
long wordendindex[wordcount];//to store wordendindex' ending indices from the mapping file
string words[wordcount]; //to store the words from mapping file.
long wordcluster[cbsize][wordcount];//to count the most number of word type vectors in a cluster, to decide which word that cluster represents

ifstream inpt,logmap;
ofstream logcep,logdist,cb, logclustersize;



struct node
{
    long cluster;
    double data[cepsize];
};

list<struct node *>nodelist,codebook; //to store the universe
list<struct node *>::iterator iter; //to iterate trhough list
list<struct node *>::iterator iternl;//to iterate through test data
list<struct node *>::iterator itercb; // to iterate through codebook vectors
struct node* centers[cbsize]; //to store sum of vectors of a cluster
//struct node* centersquare[cbsize]; //to store sum of squares of vectors of a cluster


void readceps() // this method reads all the universe of cepstral coefficients into list called nodelist
{
	string item;
    inpt.open(universe,ios::in);
    if(!inpt) //display error message if we cant open the file
    {
        cout<<"file cant be open"<<endl;
	    inpt.close();
		system("pause");
    }

	struct node* temp; //temp node to take each vector from textfile into nodelist
    while(!inpt.eof())
    {
        getline(inpt,line);//get the line by reading until newline is encountered
		cepcoecount++; //increment the number of cepstral coefficients by one
		stringstream ss(line); // convert the line to string stream so that it can be splitted easily
		int i=0;
		temp=new node;
		while(getline(ss,item,delim))//split w.r.t delim and splitted substrings lies in item
		{
			temp->data[i++]=stod(item.c_str());//c.str()-->to get a pointer to a "null-terminated character array with data equivalent to those stored in the string"
			//cout<<temp->data[i-1]<<endl;
		}
		temp->cluster=-1; //indicates not yet clustered
		nodelist.push_back(temp);
    }
	nodelist.pop_back(); //to eliminate the last line from text file which is nothing but new line and hence stores garbage results
	cepcoecount--;//so decrement the count by 1
    cout<<"number of cepstral coefficients = "<<cepcoecount<<endl;
    inpt.close();
}

void initialization() //to select initial cluster centers
{
	time_t timer;
	long ini;
	struct node* tempcopy,*tempstore;
	for(int i=0;i<cbsize;i++)
	{
		time(&timer); //get the time in ticks
		ini=(rand()+timer)%(cepcoecount-0+1)+0;//min-->0, max-->cepcoecount for random number generation between two numbers
		//cout<<ini<<endl;
		iter = nodelist.begin();
		for(long j=1;(j<ini) && (iter!=nodelist.end());j++) //skip the initial vectors to reach the desired vector
			iter++;
		tempcopy=(*iter);
		tempcopy->cluster=i; //update the cluster number
		tempstore=new node;
		tempstore->cluster=i; //assign cluster number
		for(int j=0;j<cepsize;j++) //copy the data values
			tempstore->data[j]=tempcopy->data[j];
		codebook.push_back(tempstore); //store it in codebook
		//cout<<"stored in code book "<<i<<endl;
	}
}

void initializecentroid()
{
	for(int i=0;i<cbsize;i++)
	{
		(centers[i])->cluster=0;//use cluster to store number of elements in the cluster as index already indicates the cluster number
		//(centersquare[i])->cluster=0;//use cluster to store number of elements in the cluster as index already indicates the cluster number
		for(int j=0;j<cepsize;j++)
		{
			(centers[i])->data[j]=0; //set all the data values to zero
			//(centersquare[i])->data[j]=0; //set all the data values to zero
		}
	}
}

void centroid(struct node *temp)
{
	((centers[temp->cluster])->cluster)++; //increment the number of elements in the cluster by 1
	for(int i=0;i<cepsize;i++) //add the corresponding dimensional elements which later deviding with cluster size gives centroid
		(centers[temp->cluster])->data[i]=(centers[temp->cluster])->data[i] + temp->data[i];
}

void eucledian() //classification based on eucledian distance
{
	double refdist=0,testdist=0;
	int index=0; //index of the cluster
	distortion=0;
	iternl=nodelist.begin(); //get the first vector from the universe
	while(iternl!=nodelist.end()) //iterate through universe of the vectors until all vectors are explored
	{
		itercb=codebook.begin(); //get the first code vector of the codebook
		refdist=0; //initialize reference distance
		for(int j=0;j<cepsize;j++) //assume the distance from first vector as reference distance
		{
			refdist+=((*iternl)->data[j]-(*itercb)->data[j])*((*iternl)->data[j]-(*itercb)->data[j]); //compute distance from first code vector as reference distance
		}
		index=0; //assume the vector belongs to first cluster
		for(int i=1;i<cbsize;i++)
		{
			itercb++; // go to next codebook vector
			testdist=0;
			for(int j=0;j<cepsize;j++)//find the distance from other vectors
			{
				testdist+=((*iternl)->data[j]-(*itercb)->data[j])*((*iternl)->data[j]-(*itercb)->data[j]); // compute the distance from the center of a cluster
			}
			if(testdist<refdist) //make the shortest distance as reference distance and save the correpsonding codebook vector index
			{
				refdist=testdist; //make the current lowest distance as the reference distance
				index=i; //and the vector may belongs to corresponding cluster
			}
		}
		distortion+=refdist; //add the distance to distortion
		(*iternl)->cluster=index; //associate the vector with the correpsonding cluster with lowest distance
		centroid(*iternl);//send it for centroid calculation
		//cout<<"distortion = "<<distortion<<"  cluster = "<<(*iternl)->cluster<<endl;
		iternl++; //get next vector from universe
	}
}

void tokhura()
{
	double refdist=0,testdist=0;
	int index=0; //index of the cluster
	distortion=0;
	long nodecounter=0; //to count the nodes
	memset(wordcluster, 0, sizeof wordcluster); //initialize whole 2D array to zero
	iternl=nodelist.begin(); //get the first vector from the universe
	while(iternl!=nodelist.end()) //iterate through universe of the vectors until all vectors are explored
	{
		itercb=codebook.begin(); //get the first code vector of the codebook
		refdist=0; //initialize reference distance
		for(int j=0;j<cepsize;j++) //assume the distance from first vector as reference distance
		{
			refdist+=w[j]*((*iternl)->data[j]-(*itercb)->data[j])*((*iternl)->data[j]-(*itercb)->data[j]); //compute distance from first code vector as reference distance
		}
		refdist/=cepsize;//as per tokhura distance
		index=0; //assume the vector belongs to first cluster
		for(int i=1;i<cbsize;i++)
		{
			itercb++; // go to next codebook vector
			testdist=0;
			for(int j=0;j<cepsize;j++)//find the distance from other vectors
			{
				testdist+=w[j]*((*iternl)->data[j]-(*itercb)->data[j])*((*iternl)->data[j]-(*itercb)->data[j]); // compute the distance from the center of a cluster
			}
			testdist/=cepsize; //as per tokhura distance formula
			if(testdist<refdist) //make the shortest distance as reference distance and save the correpsonding codebook vector index
			{
				refdist=testdist; //make the current lowest distance as the reference distance
				index=i; //and the vector may belongs to corresponding cluster
			}
		}
		distortion+=refdist; //add the distance to distortion
		(*iternl)->cluster=index; //associate the vector with the correpsonding cluster with lowest distance
		for(int i=0;i<wordcount;i++)
		{
			if(nodecounter<=wordendindex[i])
			{
				wordcluster[index][i]++;
				break;
			}
		}
		nodecounter++; //increment the node counter
		centroid(*iternl);//send it for centroid calculation
		//cout<<"distortion = "<<distortion<<"  cluster = "<<(*iternl)->cluster<<endl;
		iternl++; //get next vector from universe
	}
}

void tokhuraweightsbyme()
{
	struct node* temp, *xmean, *xsquare,*variance;
	temp=new node;
	xmean=new node; // to store sum of the data elements
	xsquare=new node; //to sum of the squares of the elements
	variance=new node; //to store the variance
	for(int i=0;i<cepsize;i++) //initialize the values to zero
	{
		xmean->data[i]=0;
		xsquare->data[i]=0;
	}
	iter = nodelist.begin(); //get the first vector from universe

	while(iter != nodelist.end()) //iterate till the last vector in the universe is obtained
	{
		temp=(*iter);
		for(int i=0;i<cepsize;i++) //get each data element of a vector
		{
			(xmean->data[i])+=temp->data[i]; //add the data element
			(xsquare->data[i])+=(temp->data[i])*(temp->data[i]); //add the square of the data element
		}
		iter++; //go to next vector of the universe
	}
	cout<<"tokhura weights are"<<endl;
	for(int i=0;i<cepsize;i++) //get each data element of a vector
	{
		variance->data[i]=(xsquare->data[i])/cepcoecount-((xmean->data[i])*(xmean->data[i]))/(cepcoecount*cepcoecount);
		cout << fixed << setprecision(precision); // to print with desired precision
		cout<<"SIGMA"<<i<<"="<<variance->data[i]<<"\t";
		cout<<"w["<<i<<"] = "<<1/variance->data[i]<<endl;
		w[i]=1/variance->data[i]; //reciprocal of the variance gives tokhura weight
	}
}

void classification()
{
	//eucledian(); //classification based on eucledian distance
	tokhura(); //classification based on tokhura distance
}

void updatecodevector() // to update the code vectors to centroids of the clusters
{
	long i=0;
	iterations++;
	logclustersize<<"Iteration number: "<<iterations<<endl;
	itercb=codebook.begin(); //get the first code vector of the codebook
	while(itercb!=codebook.end()) //iterate till last code vector is obtained
	{
		if(centers[i]->cluster==0) //indicates the cluster do not have any vectors associated with it
		{
			cout<<"cluster "<<i<<" is empytcell"<<endl;
			logclustersize<<"cluster "<<i<<" is empytcell"<<endl;
			itercb++; //go to next code vector
			i++;
			continue;
		}
		for(int j=0;j<cepsize;j++) //get each data element from the code vector
			(*itercb)->data[j]=(centers[i]->data[j])/(centers[i]->cluster);

		itercb++; //get next code vector
		logclustersize<<"cluster "<<i<<" ---> "<<centers[i]->cluster<<" vectors "<<endl;
		i++;
	}
	logclustersize<<endl;
	initializecentroid();//erase the old values by initializing the values to zeros
}

bool termination() //to check the termination condition
{
	if(olddistortion==0) //executed for the first iteration of the algo.
	{
		olddistortion=distortion; //assign current distortion to old distortion
		return false;
	}
	if((olddistortion-distortion)*100/distortion<thresholdist) //check the threshold condition
		return true;
	olddistortion=distortion; //store the current distortion in olddistortion for later comparison
	return false;
}

void displayvectors() //to display all vectors of cepstral coefficients
{
	struct node* temp;
	iter = nodelist.begin(); //get the first vector from universe

	while(iter != nodelist.end()) //iterate till the last vector in the universe is obtained
	{
		temp=(*iter);
		for(int i=0;i<cepsize;i++) //get each data element of a vector
			cout << temp->data[i] << " ";
		cout<<endl;
		iter++; // move to next vector
	}

}

void displaycodebook() // to display code vectors in code book
{
	struct node* temp;
	iter = codebook.begin();
	cout<<"vectors in codebook"<<endl;
	cb.open("codebook.txt");
	if(!cb) //if we cant open the file notify with error message
	{
		cout<<"can't open codebook.txt file"<<endl;
		cb.close();
		return;
	}
	while(iter != codebook.end()) //iterate till the end of the codebook
	{
		temp=(*iter);
		for(int i=0;i<cepsize;i++) //get each data element from the vector
		{
			cout << fixed << setprecision(precision); // to print with desired precision
			cout << temp->data[i] << " "; //write to standard output
			cb << fixed << setprecision(precision); //to store with desired precision
			cb<<temp->data[i] << " "; //write to file
		}
		cout<<endl;
		cb<<endl;
		iter++; //move to next code vector
	}
	cb.close(); //close the file
}

void disp_cluster_repres() //display which word each cluster represents
{
	cout << "\n\n    cluster\tword\t" <<endl;
	for(int i=0;i<cbsize;i++)
	{
		int max=0;
		for(int j=0;j<wordcount;j++)
		{
			if(wordcluster[i][max]<wordcluster[i][j])
				max=j;
		}
		if(wordcluster[i][max]==0)
		{
			cout<<"\t"<<i << "\tempty cell"<<endl;
				continue;
		}
		cout <<"\t"<< i<<"\t" <<words[max]<<endl;
	}
}

void tokhuraweightsbysir()
{
	w[0]=1;//these are tokhura weights given by sir
	w[1]=3;
	w[2]=7;
	w[3]=13;
	w[4]=19;
	w[5]=22;
	w[6]=25;
	w[7]=33;
	w[8]=42;
	w[9]=50;
	w[10]=56;
	w[11]=61;
}

void algo()
{
	for(int i=0;i<cbsize;i++) //create it once and initialize after every iteration
	{
		centers[i]=new node;
		//centersquare[i]=new node;
	}
	initialization(); //to select initial cluster centers
	initializecentroid(); // to initialize all centroid values to zeros

	logdist.open("logdistortion.txt"); //log file to store the distortion after each iteration
	logclustersize.open("logclustersize.txt");
	while(!termination())
	{
		classification(); // classify the vectors based on nearest cluster center
		updatecodevector(); //update the cluster centers to new centroids
		cout<<"distortion = "<<distortion<<"\titeration = "<<iterations<<endl;
		logdist<<distortion/cbsize<<endl; //average distortion per codebook vector or per cluster
		//logdist<<abs((olddistortion-distortion))*100/distortion<<endl;
	}
	logdist.close();
	logclustersize.close();
}

void kmeans()
{
	readceps(); //to read the cepstral coefficient vectors from available universe
	tokhuraweightsbysir(); //use the tokhura weights given by sir
	//tokhuraweightsbyme(); //compute the tokhura weights from universe
	algo(); // apply algorithm on universe
	//displayvectors(); //diplay vectors of universe
	//displaycodebook(); //store and display generated codebook
}

int main()
{
	string splitmap;//used to split each line in mapping file with respect to mapdelim delimiter
	logmap.open("mapping.txt",ios::in);//read the mapping from vowels to vectors
	if(!logmap) //display error message if we cant open the file
    {
        cout<<"file can not be open"<<endl;
	    logmap.close();
		system("pause");
		return 0;
    }
	getline(logmap, line); //read and ignore first line
	int temp=0; //to count the lines in mapping file
	cout<<"reading mapping file..."<<endl;
	while(!logmap.eof())
	{
		string tempstring; //to read temporary string which will further be converted into long
		getline(logmap,line);
		stringstream ss(line);
		int j=0;
		while(getline(ss,splitmap,mapdelim))//split w.r.t delim and splitted substrings lies in splitmap
		{
			if(j==0)
				words[temp++]=splitmap;
			if(j!=2) //ending index of each word is stored in second column of mapping file.
			{
				j++;
				continue;
			}
			j++;
			wordendindex[temp]=atol(splitmap.c_str());
			//cout<<splitmap->data[i-1]<<endl;
			cout<<temp-1<<"th word ending index is "<<wordendindex[temp-1]<<endl;
		}
	}
	logmap.close();
	kmeans(); //call k-means
	disp_cluster_repres(); //display which word each cluster represents
	system("pause");
	//delete []wordcluster;
	return 0;
}
