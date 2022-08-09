#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string.h>
#include <sstream>
#include <list>
#include <conio.h>
#include "config.h"

using namespace std;

double pi[N];
double a[N][N];
double b[N][M];
double pi_bar[N];
double a_bar[N][N];
double b_bar[N][M];
double P; //to store the probability in apha, beta, gamma computations
double **alpha;
double **beta;
double **gamma;
double **delta;
int **psi;
double ***zeta;
double Pstar; //maximum probability in viterbi algo.
int *Qstar;//

std::ostringstream oss;

string filename, line;
short samples[600000], MaxAmp=0; //Maximum number of samples we can handle and Maximum amplitude desired for scaling
double ThresholdZCR, DCshift=0; //Initialize DCShift and Threshold ZCR
double TotalEnergy=0, ThresholdEnergy, NoiseEnergy=0, TotalZCR=0; //TotalZCR is only for first few InitFrames
long start=0, stop=0, framecount, samplecount=0; //start and end marker for speech signal
string filenam;
unsigned short q; //number of coefficients (ai's) that need to be found
long double *R, *c, *E, *k, *alfa;
long double *x;
double w[cepsize];//to store tokhura weights given by sir

int T; //to store number of time frames
int *o; //to store the observation sequence

ifstream InpSpeech,cb;
ofstream ScaledInpSpeech,out,cep,logcep,logmap;

struct node //to represent each vector
{
    long cluster;
    double data[cepsize];
};
list<struct node *>cblist;
list<struct node *>testlist; //to store cepstral coefficients of test speech
list<struct node *>::iterator itertest, itercb; //to iterate trhough list

void durbin() //apply durbin method in order to compute the coefficients
{
	E = new long double[p+1];
	k = new long double[p+1];
	alfa = new long double[p+1]; // for new alpha values
	long double* alp = new long double[p+1]; // for old values
	if(R[0]==0)
	{
		cout << "\nEnergy should not be ZERO\n";
		return;
	}
	E[0]=R[0]; //Initialize E[0] with energy
	for(int i=0;i<=p;i++)
	{
		alfa[i]=0; //initialize all the initial coefficients to 0
		alp[i]=0;
	}
	for(int i=1;i<=p;i++)
	{
		if(i==1)
			k[1]=R[1]/R[0]; //special case when i=1
		else //find k(i) for all other remaining values
		{
			long double sum=0;
			for(int j=1;j<=i-1;j++)
			{
				alp[j]=alfa[j]; // store the old values in alp array
				sum+=alp[j]*R[i-j];
			}
			k[i]=(R[i]-sum)/E[i-1];
		}
		alfa[i]=k[i]; //assign values to coefficients
		for(int j=1;j<=i-1;j++)
			alfa[j]=alp[j]-k[i]*alp[i-j]; //update coefficients from previously obtained values
		E[i]=(1-k[i]*k[i])*E[i-1]; //update E(i)
	}
/*	cout << "\nLPC Coefficient values\n";
	for(int i=1;i<=p;i++)
		cout << "a" << i << " = " << alfa[i] << "\n";
*/
}

void cepstralcoefficients() //calculate the cepstral coefficients
{
	unsigned int m,j;
	c = new long double[q+1];
	struct node* temp;
	temp=new node;

	c[0]= logl(R[0]); //initial cepstral coefficent computation from energy
	for(m=1;m<=p;m++)
	{
		long double sum=0;
		for(j=1;j<=m-1;j++) //calculate the sum from older cepstral coefficents to compute new cepstral coefficients
			sum+=j*c[j]*alfa[m-j]/m;
		c[m]=alfa[m]+sum; //new cepstral coefficients
	}
	if(m>p) // for our assignment this never get executed as we assume q=p
	{
		for(;m<=q;m++)
		{
			long double sum=0;
			for(j=m-p;j<=m-1;j++)
				sum+=j*c[j]*alfa[m-j]/m;
			c[m]=sum;
		}
	}
	for(int i=1;i<=q;i++)
	{   
		temp->data[i-1]= c[i]; //we do not store c[0], and when we read back the codebook c[1] comes to data[0]
	}
	testlist.push_back(temp);
}

void Coefficients()
{

	q=p; // for our assignment we always take q=p

	durbin(); //call this function to evaluate durbin algorithm
	cepstralcoefficients(); //call method to calculate cepstral coefficients

}

void AutoCorrelation() //calculate the auto correlation
{
	long n=framesize, ind=0;
	R = new long double[p+1]; //create array of size p+1 for p coefficients
	for(ind=0;ind<=p;ind++)  //generate p+1 R values to compute cepstral coefficients
	{
		double square=0;
		for(int t=0;t<framesize-ind;t++)
		{
			square+=x[t]*x[t+ind];
		}
		R[ind]=square/n;
	}
	//cout << "frame number " << K <<endl;
	Coefficients(); //call this function to compute the coefficients
}

bool Repeat(int slipcount) //repeated for each frame to calculate the cepstral coefficients
{
		samplecount=0;
		x=new long double[framesize]; // to store all the samples of a frame
		InpSpeech.open(testt, ios::in); // open the file to read from it
		if (!InpSpeech) //if file is not present pop the error
        {
    		cout << "\n**File failed to open**\n\n";
            InpSpeech.clear();
        }

		//count the number of samples and frames in the file
		if (InpSpeech.is_open()) //if file is open then only execute this block
		{
			while ( !InpSpeech.eof() ) //read till the end of file is reached
			{
				getline (InpSpeech,line); //read a line at a time
				samplecount+=1;     //increment the sample count by 1
				if(samplecount > IgnoreSamples*slipcount) //first 5 lines in text file indiactes type of encoding of the speech,
				{
					x[samplecount - (IgnoreSamples*slipcount + 1)] = (double)atof(line.c_str()); 
					if(framesize==samplecount - (IgnoreSamples*slipcount + 1))
						break;
				}
			}
			InpSpeech.close();
		}
		samplecount = samplecount - (IgnoreSamples*slipcount + 1);
		if(samplecount < framesize)
		{
			logcep << "less number of samples "<<samplecount<<endl<<endl;
			return false;
		}

		AutoCorrelation();
		//Coefficients(); //call this function to compute the coefficients

		//cout <<"open cep.txt to see the cepstral coefficients for each frame" <<endl;
		return true;
}

void readcepsfromcb() // this method reads all the codebook of cepstral coefficients into list called cblist
{
	string item;
    cb.open(codebook,ios::in);
    if(!cb) //display error message if we cant open the file
    {
        cout<<"file cant be open"<<endl;
	    cb.close();
		system("pause");
    }

	struct node* temp; //temp node to take each vector from textfile into cblist
    while(!cb.eof())
    {
        getline(cb,line);//get the line by reading until newline is encountered
		stringstream ss(line); // convert the line to string stream so that it can be splitted easily
		int i=0;
		temp=new node;
		while(getline(ss,item,delim))//split w.r.t delim and splitted substrings lies in item
		{
			temp->data[i++]=stod(item.c_str());//c.str()-->to get a pointer to a "null-terminated character array with data equivalent to those stored in the string"
			//cout<<temp->data[i-1]<<endl;
		}
		temp->cluster=-1; //indicates not yet clustered
		cblist.push_back(temp);
    }
	//cblist.pop_back(); //to eliminate the last line from text file which is nothing but new line and hence stores garbage results
    cb.close();
}

void tokhuraweightsbysir() //these are tokhura weights given by sir
{
	w[0]=1;
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

void observation_sequence() //find the observation sequence based on shortest distance
{
	double refdist=0,testdist=0;
	int index=0,tempcount=0;
	tokhuraweightsbysir();
	itertest=testlist.begin();
	while(itertest!=testlist.end())
	{
		itercb=cblist.begin();
		refdist=0; //initialize reference distance
		for(int j=0;j<cepsize;j++) //assume the distance from first vector as reference distance
		{
			refdist+=w[j]*((*itertest)->data[j]-(*itercb)->data[j])*((*itertest)->data[j]-(*itercb)->data[j]); //compute distance from first code vector as reference distance
		}
		refdist/=cepsize;//as per tokhura distance
		index=0; //assume the vector belongs to first cluster
		for(int i=1;i<cbsize;i++)
		{
			itercb++; // go to next codebook vector
			testdist=0;
			for(int j=0;j<cepsize;j++)//find the distance from other vectors
			{
				testdist+=w[j]*((*itertest)->data[j]-(*itercb)->data[j])*((*itertest)->data[j]-(*itercb)->data[j]); // compute the distance from the center of a cluster
			}
			testdist/=cepsize; //as per tokhura distance formula
			if(testdist<refdist) //make the shortest distance as reference distance and save the correpsonding codebook vector index
			{
				refdist=testdist; //make the current lowest distance as the reference distance
				index=i; //and the vector may belongs to corresponding cluster
			}
		}
		o[tempcount]=index;
		itertest++;
		tempcount++;
	}
}

void disp_obs_seq() //display observation sequence
{
	cout<< "\n================observation sequence============="<<endl;
	out<< "\n================observation sequence============="<<endl;
	for(int i=0;i<T;i++)
	{
		cout<<o[i]<<"\n";
		out<<o[i]<<"\n";
	}
}

void diplay_model(double pi[N], double a[N][N], double b[N][M]) //display the whole model
{
	cout<<"Initial State probabilities: "<<endl;
	out<<"Initial State probabilities: "<<endl;
	for(int i=0;i<N;i++)  //display initial state state distributions
	{
        cout<<pi[i]<<"\t";
		out<<pi[i]<<"\t";
	}
	cout<<"\ntransition probabilities: "<<endl;
	out<<"\ntransition probabilities: "<<endl;
    for(int i=0;i<N;i++) //display state transition probabilities
	{
        for(int j=0;j<N;j++)
		{
            cout<<a[i][j]<<"\t";
			out<<a[i][j]<<"\t";
		}
		cout<<endl;
		out<<endl;
	}
	cout<<"observation symbol probability: "<<endl;
	out<<"observation symbol probability: "<<endl;
    for(int i=0;i<N;i++) //display observation symbol probability distribution
	{
        for(int j=0;j<M;j++)
		{
            cout<<b[i][j]<<"\t";
			out<<b[i][j]<<"\t";
		}
		cout<<endl;
		out<<endl;
	}
}

void uniform_model() //all transition probabilities and observation symbol probabilities are uniformly distributed
{
	//cout<<"In forward procedure"<<endl;
	//out<<"In forward procedure"<<endl;
	for(int i=0;i<N;i++) //assign the given values
        pi[i]=1.0/N;
    for(int i=0;i<N;i++) //assign the given values for transition probability distribution
        for(int j=0;j<N;j++)
            a[i][j]=1.0/N;
    for(int i=0;i<N;i++) //assign the given values for observation symbol probability distribution
        for(int j=0;j<M;j++)
            b[i][j]=1.0/M;
}

void feed_forward_model() //this is biased model which works in lazy way
{
	//cout<<"In forward procedure"<<endl;
	//out<<"In forward procedure"<<endl;
	for(int i=0;i<N;i++) //assign the given values
		if(i==0) //make the first state as starting state
			pi[i]=1.0;
		else
			pi[i]=0;
    for(int i=0;i<N;i++) //assign the given values for transition probability distribution
        for(int j=0;j<N;j++)
			if(i==j&&i!=N-1)
				a[i][j]=0.8; //probability of being in current state
			else if(i==j&&i==N-1)
				a[i][j]=1; //forcing to end the transition in final state
			else if(j==i+1)
				a[i][j]=0.2; //probability to transit to next immediate state
			else
				a[i][j]=0; //probability to move to non immediate states
    for(int i=0;i<N;i++) //assign the given values for observation symbol probability distribution
        for(int j=0;j<M;j++)
            b[i][j]=1.0/M;
}

double forward_proc() //alpha computation
{
    for(int i=0;i<N;i++) //initialization
        alpha[0][i]=pi[i]*b[i][o[0]];
    for(int t=0;t<T-1;t++) //induction
    {
        for(int j=0;j<N;j++)
        {
            double sum=0;
            for(int i=0;i<N;i++)
            {
                sum+=alpha[t][i]*a[i][j];
            }
            alpha[t+1][j]=sum*b[j][o[t+1]];
        }
    }
    P=0;
    for(int i=0;i<N;i++) //estimate what is the probability that the observation sequence is from the current model
	{
        P+=alpha[T-1][i];
		//cout<<alpha[T-1][i]<<endl;
	}
    return P;
}

double backward_proc() //beta computation
{
    for(int i=0;i<N;i++) //initialization
        beta[T-1][i]=1;
    for(int t=T-2;t>=0;t--) //induction
    {
        for(int i=0;i<N;i++)
        {
            double sum=0;
            for(int j=0;j<N;j++)
            {
                sum+=a[i][j]*b[j][o[t+1]]*beta[t+1][j];
            }
            beta[t][i]=sum;
        }
    }
    P=0;
    for(int i=0;i<N;i++) //this is for self assessment and not required as there is no further use
    {
        P+=beta[0][i];
		//cout<<beta[0][i]<<endl;
    }
   return P;
}

double gamma_proc() //gamma computation
{
	int *q,argmax=0; //q--> store the state which has maximum probability of occurence at time t.
	q=new int[T];
	double devider=0; //used as devider in baye's theorem for computation of gamma
	for(int t=0;t<T;t++)
	{
		for(int i=0;i<N;i++) //compute it once for t
		{
			devider+=alpha[t][i]*beta[t][i];
		}
		argmax=0;
		for(int i=0;i<N;i++)
		{
			gamma[t][i]=alpha[t][i]*beta[t][i]/devider;
			if(gamma[t][argmax]<gamma[t][i])
				argmax=i;
		}
		q[t]=argmax;
		devider=0;
	}
	P=1;
	for(int t=0;t<T;t++)
	{
		P*=gamma[t][q[t]];
		//cout<<gamma[t][q[t]]<<endl;
	}
	return P;
}

double viterbi_proc() //delta and psi computation
{
	Qstar=new int[T];
	int argmax=0; // to store the argument which gives maximum probability
	for(int i=0;i<N;i++) //initialization
	{
		delta[0][i]=pi[i]*b[i][o[0]];
		//cout <<delta[0][i]<<endl;
		psi[0][i]=-1; //indicates, no state is yet assigned
	}
	for(int t=1;t<T;t++) //recursion over time sequence
	{
		for(int j=0;j<N;j++) // to store maximum probabilities
		{
			argmax=0; //assume the first state gives maximum probability
			for(int i=1;i<N;i++)
			{
				if(delta[t-1][i]*a[i][j] > delta[t-1][argmax]*a[argmax][j]) //checking for maximum score
					argmax=i; //argument which gives maximum score
			}
			delta[t][j]=delta[t-1][argmax]*a[argmax][j]*b[j][o[t]];
			psi[t][j]=argmax; //largest argument index which gives maximum probability
		}
	}
	argmax=0;
	for(int i=1;i<N;i++) //to find the argmax for last time frame
		if(delta[T-1][i] > delta[T-1][argmax])
			argmax=i;
	Pstar=delta[T-1][argmax];
	//cout<<"state sequence with back tracking"<<endl;
	Qstar[T-1]=argmax;
	//cout<<T-1+1<<" ---> "<<Qstar[T-1]+1<<endl;
	for(int t=T-2;t>=0;t--) //back tracking the path
	{
		Qstar[t]=psi[t+1][Qstar[t+1]];
		//cout<<t+1<<" --> "<<Qstar[t]+1<<endl;
	}
	return Pstar;
}

void disp_state_seq() //display the state sequence followed
{
	cout<<"\nobserved state sequence: "<<endl;
	out<<"\nobserved state sequence: "<<endl;
	for(int t=0;t<T;t++)
	{
		cout<<Qstar[t]<<endl;
		out<<Qstar[t]<<endl;
	}
}

void baum_welch_proc() //zeta computaion
{
	double devider=0; //used as common devider just like in baye's theorem
	for(int t=0;t<T-1;t++) //repeat this for all T-1 state transitions
	{
		devider=0;
		for(int i=0;i<N;i++)
		{
			for(int j=0;j<N;j++)
				devider+=alpha[t][i]*a[i][j]*b[j][o[t+1]]*beta[t+1][j];
		}
		for(int i=0;i<N;i++)
		{
			for(int j=0;j<N;j++)
				zeta[t][i][j]=alpha[t][i]*a[i][j]*b[j][o[t+1]]*beta[t+1][j]/devider;
		}
	}
}

void re_estimation() //re-estimate transition probabilities and observation symbol probabilities
{
	double numerator, denominator; //for re-estimation of transition probabilities
	for(int i=0;i<N;i++) //re-estimation of pi as pi_bar
		pi_bar[i]=gamma[0][i];
	for(int i=0;i<N;i++) //re-estimation of a as a_bar
	{
		for(int j=0;j<N;j++)
		{
			numerator=0;
			denominator=0;
			for(int t=0;t<T-2;t++)
			{
				numerator+=zeta[t][i][j];
				denominator+=gamma[t][i];
			}
			a_bar[i][j]=numerator/denominator;
		}
	}
	for(int j=0;j<N;j++) //re-estimation of b as b_bar
	{
		for(int k=0;k<M;k++)
		{
			numerator=0;
			denominator=0;
			for(int t=0;t<T;t++)
			{
				if(o[t]==k)
					numerator+=gamma[t][j];
			}
			for(int t=0;t<T-1;t++)
			{
				denominator+=gamma[t][j];
			}
			b_bar[j][k]=numerator/denominator;
		}
	}
}

int read_test() //read the test input from the speaker through mic(if desired)
{
	long i,j;
	cout<<"Do you want to record a speech?(y/n): ";
	switch(_getch())
	{
		case 'y':
					oss << recmod << " "<<duration <<" "<<testw << " " << testt; //define the string manually, because macro cant expand inside string
					system(oss.str().c_str()); //call the system function to read directly from the mic
					cout<<"\nyour speech has been succesfully recorded"<<endl;
					break;
		case 'n':	cout<< "\nwe will try to recognize the speech in "<<testt<<" file"<<endl; //execute this condition when the offline test sspeech is already available
					break;
		default:	cout<< "\nSorry, no such option is entertained"<<endl;
					cout<<"BYE.... BYE...."<<endl;
					system("pause");
					return 0;
	}
	filename=testt; //name of the recorded speech text file
	InpSpeech.open(filename, ios::in); // open the file to read from it
	if (!InpSpeech) //if file is not present pop the error
    {
   		cout << "\n**File failed to open**\n\n";
        InpSpeech.clear();
    }
	out.open("out.txt"); // open the file to write the results
	//count the number of samples and frames in the file
	if (InpSpeech.is_open()) //if file is open then only execute this block
	{
		while ( !InpSpeech.eof() ) //read till the end of file is reached
		{ 
			getline (InpSpeech,line); //read a line at a time
			samplecount+=1;     //increment the sample count by 1
			if(samplecount > IgnoreSamples + 4) //first 4 lines in text file indiactes type of encoding of the speech,
			{
				samples[samplecount - (IgnoreSamples + 5)] = (short)atoi(line.c_str()); //5=4+1 4-->indicates encoding, 1-->array index starts with 0
				DCshift+=samples[samplecount - (IgnoreSamples + 5)];
				if(abs(samples[samplecount - (IgnoreSamples + 5)])>MaxAmp)
					MaxAmp=abs(samples[samplecount - (IgnoreSamples + 5)]); //MaxAmp contains the magnitude od the maximum valued(absolute) sample
			}
		}
		InpSpeech.close();
	}
	samplecount = samplecount - (IgnoreSamples + 4);
	framecount = samplecount/framesize;
	DCshift=DCshift/samplecount;
	out << "Number of Samples = " << samplecount << "\n"; // print the results in out.txt file
	out << "Number of frames = " << framecount <<" each of size "<<framesize<< " samples\n";
	out << "DC shift needed is " << DCshift << "\n";
	out << "Maximum Amplitude = " << MaxAmp << "\n";
	cout << "Number of Samples = " << samplecount << "\n"; // print the results on standard output
	cout << "Number of frames = " << framecount <<" each of size "<<framesize<< " samples\n";
	cout << "DC shift needed is " << DCshift << "\n";
	cout << "Maximum Amplitude = " << MaxAmp << "\n";
	
	cout<<"Did you remove the non speech part of the signal?(y/n): ";
	switch(_getch())
	{
		case 'n':
					cout<<"\nNo problem, we will remove it for you..." <<endl;
					//Store scaled samples in different file
					ScaledInpSpeech.open("ScaledInpSpeech.txt"); // open ScaledInpSpeech.txt to write the scaled sample values
					for(i=0;i<samplecount;i++)
						ScaledInpSpeech << (samples[i] - DCshift)*Ampscale/MaxAmp << "\n"; // writing the scaled samples to the file ScaledInpSpeech.txt
					ScaledInpSpeech.close();
					//use the scaled samples
					InpSpeech.open("ScaledInpSpeech.txt", ios::in); // open ScaledInpSpeech.txt to read the scaled samples
					if (!InpSpeech)
					{
						cout << "\n**File failed to open**\n\n";
					    InpSpeech.clear();
					}
					InpSpeech.close(); // close the file
					//ZCR and Energy calculation
					double AvgZCR[MaxFrameCount], AvgEnergy[MaxFrameCount];
					for(i=0;i<framecount;i++)
					{
						AvgZCR[i]=0; // Initialize Average ZCR for each frame
						AvgEnergy[i]=0; // Initialize Average Energy for each frame
						for(j=0;j<framesize-1;j++)
						{
							if((samples[i*framesize + j]-DCshift)*Ampscale/MaxAmp*((samples[i*framesize + j + 1]-DCshift)*Ampscale*1.0 < 0)) // If two adjacent samples have opposite sign then increment the avergae ZCR by 1
								AvgZCR[i]+=1;
							AvgEnergy[i]+=1.0*(samples[i*framesize + j]-DCshift)*Ampscale/MaxAmp*(samples[i*framesize + j]-DCshift)*Ampscale/MaxAmp; // energy is calculated as square of amplitude of the sample value
						}
						//out << "ZCR in "<<i+1<<"th frame is "<< AvgZCR[i] << "\n";
						AvgZCR[i]/=framesize;
						//out << "Avergae ZCR in\t\t "<<i+1<<"th frame is "<< AvgZCR[i] << "\n";
						AvgEnergy[i]+=1.0*(samples[i*framesize + j]-DCshift)*Ampscale/MaxAmp*(samples[i*framesize + j]-DCshift)*Ampscale/MaxAmp; // average energy for all the frames after scaled sample values
						//out << "Energy in\t "<<i+1<<"th frame is "<< AvgEnergy[i] << "\n";
						AvgEnergy[i]/=framesize; // average energy per frame after scaled sample values 
						//out << "Average Energy in \t\t\t"<<i+1<<"th frame is "<< AvgEnergy[i] << "\n";
						TotalEnergy+=AvgEnergy[i];
					}
					//calculate noise energy to decide threshold for energy
					for(i=0;i<InitFrames;i++)
					{
						TotalZCR+=AvgZCR[i];
						NoiseEnergy+=AvgEnergy[i];
					}
					ThresholdZCR=TotalZCR/InitFrames;
					NoiseEnergy/=InitFrames;
					ThresholdZCR*=0.9;
					//ThresholdEnergy=TotalEnergy/framecount; //threshold energy is the average of all averaged energies
					//ThresholdEnergy*=0.9;
					ThresholdEnergy=NoiseEnergy*10;
					bool flag;
					flag=false;
					//start and end marker of speech
					for(i=0;i<framecount-3;i++)
					{
						//if(AvgEnergy[i+1]>ThresholdEnergy)
						if(AvgZCR[i+1]<ThresholdZCR || AvgEnergy[i+1]>ThresholdEnergy)
						{
							if(flag == false && AvgEnergy[i+2]>ThresholdEnergy && AvgEnergy[i+3]>ThresholdEnergy)
							{
								start = i  ; //i th frame is the starting frame marker
								out << "Starting frame is "<< start+1 <<"th frame and starting sample is "<< (start+1)*framesize <<"\n"; // write starting frame and sample number in out.txt file
								out << "Starting time = " << 1.0*(start + 1)*framesize/samplingrate << " seconds\n" ;// writing starting sample time in seconds to out.txt
								cout << "Starting frame is "<< start+1 <<"th frame and starting sample is "<< (start+1)*framesize <<"\n"; // write starting frame and sample number on standard output
								cout << "Starting time = " << 1.0*(start + 1)*framesize/samplingrate << " seconds\n" ; // writing starting sample time in seconds to standard output
								flag = true;
							}
						}
						else if(flag == true && AvgZCR[i] > ThresholdZCR && AvgEnergy[i] < ThresholdEnergy && AvgEnergy[i-1] < ThresholdEnergy && AvgEnergy[i-2] < ThresholdEnergy)
						{
								stop = i ;
								out << "Ending frame is "<< stop+1 <<"th frame and Ending sample is "<< (stop+1)*framesize<<"\n"; // write ending frame and sample number in out.txt file
								out << "Ending time = " << 1.0*(stop + 1)*framesize/samplingrate << " seconds\n" ; // writing ending sample time in seconds to out.txt
								cout << "Ending frame is "<< stop+1 <<"th frame and Ending sample is "<< (stop+1)*framesize<<"\n"; // write ending frame and sample number on standard output
								cout << "Ending time = " << 1.0*(stop + 1)*framesize/samplingrate << " seconds\n" ; // writing ending sample time in seconds to standard output
								flag = false;
								break;
						}
					}
					ScaledInpSpeech.open(testt); //write the speech samples back to input file
					for(long i=start*framesize;i<=stop*framesize;i++)
					{
						ScaledInpSpeech<<(samples[i]-DCshift)*Ampscale/MaxAmp<<endl;
					}
					ScaledInpSpeech.close();
					cout << "NOTE: To see the speech samples, check ScaledInpSpeech.txt file \n";
					cout << "Check the output in out.txt file\n"; 
					break;
		case 'y': 
					start=0;
					stop= framecount;
					break;
		default:	cout<<"\nSorry, no such option vailable\n We are aborting..." <<endl;
					cout<<"BYE.... BYE...."<<endl;
					system("pause");
					return 0;
	}
	out.close();
	return 0;
}

void replace_model() //replace old model (lambda) with new model (lambda-dash)
{
	cout<<"replacing the old model with new model"<<endl;
	out<<"replacing the old model with new model"<<endl;
	for(int i=0;i<N;i++) //assign the given values
        pi[i]=pi_bar[i];
    for(int i=0;i<N;i++) //assign the given values for transition probability distribution
        for(int j=0;j<N;j++)
            a[i][j]=a_bar[i][j];
    for(int i=0;i<N;i++) //assign the given values for observation symbol probability distribution
        for(int j=0;j<M;j++)
            b[i][j]=b_bar[i][j];
}

int main()
{
	///////////////////////////////////////////////////////to read the test input signal////////////////////////////////////////////////////////////////
	read_test();

	///////////////////////////////////////////////////////////to compute the cepstral coefficients and generate codebook vector indices//////////////////////////////////////
	out.open("log.txt");
	readcepsfromcb(); //to read codebook vectors
	int skipcount=0;
	while(Repeat(skipcount++));
	T=skipcount-1;
	o = new int[T];
											//////////////////All declarations///////////////////////////////////////
	alpha=new double*[T];
	for(int i=0;i<T;i++)
		alpha[i]=new double[N];
	beta=new double*[T];
	for(int i=0;i<T;i++)
		beta[i]=new double[N];
	gamma=new double*[T];
	for(int i=0;i<T;i++)
		gamma[i]=new double[N];
	delta=new double*[T];
	for(int i=0;i<T;i++)
		delta[i]=new double[N];
	psi=new int*[T];
	for(int i=0;i<T;i++)
		psi[i]=new int[N];
	zeta=new double**[T];
	for(int i=0;i<T;i++)
	{
		zeta[i]=new double*[N];
		for(int j=0;j<N;j++)
			zeta[i][j]=new double[N];
	}
	cout<<"\nNumber of frames "<<skipcount-1<<" with skipcount = "<<IgnoreSamples<<endl;
	out<<"\nNumber of frames "<<skipcount-1<<" with skipcount = "<<IgnoreSamples<<endl;
	observation_sequence();
	disp_obs_seq();
	//uniform_model(); //all are equally probable
	feed_forward_model(); //biased transitions

	///////////////////////////////////////////////////////////alpha, beta, gamma, zeta, psi computation//////////////////////////////////////////////////////////
	int counter=0;
	while(counter++<ModelIterate)
	{
		cout<<"Iteration = "<<counter<<endl;
		cout << "Probability of an utterance being from the model: "<<forward_proc()<<endl;
		out << "Probability of an utterance being from the model: "<<forward_proc()<<endl;
		//cout << "Probability in backward procedure: "<<backward_proc()<<endl;
		//out << "Probability in backward procedure: "<<backward_proc()<<endl;
		backward_proc();
		//cout << "Probability in gamma procedure: "<<gamma_proc()<<endl;
		//out << "Probability in gamma procedure: "<<gamma_proc()<<endl;
		gamma_proc();
		cout << "Probability in viterbi procedure: "<<viterbi_proc()<<endl;
		out << "Probability in viterbi procedure: "<<viterbi_proc()<<endl;
		disp_state_seq();

		/////////////////////////////////////////////////////re-estimation solution///////////////////////////////////////////////////////////
		baum_welch_proc();
		re_estimation();

		cout<<"\n=====================MODEL before re-estimation================================"<<endl;
		out<<"\n=====================MODEL before re-estimation================================"<<endl;
		diplay_model(pi, a, b);
		cout<<"\n=====================MODEL after re-estimation================================"<<endl;
		out<<"\n=====================MODEL after re-estimation================================"<<endl;
		diplay_model(pi_bar, a_bar, b_bar);
		replace_model();
	}
	out.close();
	out.open("output.txt"); //store the final results in output.txt file
	out<<"\n=====================FINAL MODEL after re-estimation================================"<<endl;
	diplay_model(pi_bar, a_bar, b_bar);
	out<<"======================================================================================";
	out<<"\n\n";
	out << "Probability of an utterance being from the model: "<<forward_proc()<<endl;
	backward_proc();
	gamma_proc();
	disp_state_seq();
	out.close();
	system("pause");
	return 0;
}