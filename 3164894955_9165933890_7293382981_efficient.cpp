#include<bits/stdc++.h>
#include<fstream>
#include<chrono>
#include <sys/resource.h>

#define   RUSAGE_SELF     0
#define   RUSAGE_CHILDREN     -1

using namespace std;
using namespace std::chrono;
int charValue(char c)
{
	if(c=='A') return 0;
	else if (c=='C') return 1;
	else if (c=='G') return 2;
	return 3;
}


pair<int, pair<string, string> > findOptimalCost(string s1, string s2, int delta, int alphas[4][4])
{
	int m = s1.length();
	int n = s2.length();
	vector< vector<pair<int, pair<int, int> > > > dp(m+1, vector<pair<int, pair<int, int> > >(n+1));
	for(int i=0;i<=m;i++)
	{
		vector<pair<int, pair<int, int> > > vec;
		for (int j= 0; j <=n; ++j)
		{
			vec.push_back(make_pair(0, make_pair(-1,-1)));
		}
		dp.push_back(vec);
	}

	for(int i=0;i<=m;i++)
	{
		dp[i][0].first = (i)*delta;
	}

	for(int i=0;i<=n;i++)
	{
		dp[0][i].first = (i)*delta;
	}

	for(int i=1;i<=m;i++)
	{
		for(int j=1;j<=n;j++)
		{
			 dp[i][j].first = dp[i-1][j-1].first + alphas[charValue(s1[i-1])][charValue(s2[j-1])];
			 dp[i][j].second = make_pair(i-1, j-1);

			 if(dp[i-1][j].first + delta < dp[i][j].first)
			 {
			 	dp[i][j].first = dp[i-1][j].first + delta;
			 	dp[i][j].second = make_pair(i-1, j);
			 }

			 if(dp[i][j-1].first + delta < dp[i][j].first)
			 {
			 	dp[i][j].first = dp[i][j-1].first + delta;
			 	dp[i][j].second = make_pair(i, j-1);
			 }
		}
	}

	int minCostForAlignment = dp[m][n].first;

	string xAlignment = "";
	string yAlignment = "";
	int i=m,j=n;
	while(i>0 && j > 0)
	{
		if(dp[i][j].second.first == i-1 && dp[i][j].second.second == j-1)
		{
			xAlignment = s1[i-1] + xAlignment;
			yAlignment = s2[j-1] + yAlignment;
			i--;
			j--;
		}
		else if(dp[i][j].second.first == i-1 )
		{
			xAlignment = s1[i-1] + xAlignment;
			yAlignment = "_" + yAlignment;
			i--;
		}
		else
		{
			xAlignment = "_" + xAlignment;
			yAlignment = s2[j-1] + yAlignment;
			j--;
		}
	}
	while(i>0)
	{
		xAlignment = s1[i-1] + xAlignment;
		yAlignment = "_" + yAlignment;
		i--;
	}

	while(j>0)
	{
		xAlignment = "_" + xAlignment;
		yAlignment = s2[j-1] + yAlignment;
		j--;
	}

	return make_pair(minCostForAlignment, make_pair(xAlignment, yAlignment));
}

// void printVector(vector<int> vec)
// {
// 	for(int i=0;i<vec.size();i++)
// 	{
// 		cout<<vec[i]<<" ";
// 	}
// 	cout<<endl;
// 	return;
// }


vector<int> optimalCostValueOptimisedSpace(string s1, string s2, int delta, int alphas[4][4])
{
	int m = s1.length();
	int n = s2.length();

	if(n==0) 
	{
		vector<int> col1;
		col1.push_back(m*delta);
		return col1;
		// return delta * m;
	}

	if(m==0){
		vector<int> col1;
		for(int i=0;i<=n;i++)
			col1.push_back(i*delta);
		return col1;
	 	// return delta * n;
	}


	vector<int> col1(m+1);
	vector<int> col2(m+1);

	for(int i=0;i<=m;i++)
	{
		col1[i] = (i*delta);
	}
	// printVector(col1);
	vector<int> ans;
	ans.push_back(col1[m]);

	for(int i=1;i<=n;i++)
	{
		col2[0]=(i*delta);
		for(int j=1;j<=m;j++)
		{
			col2[j]= min(col2[j-1] + delta, min(col1[j-1] + alphas[charValue(s1[j-1])][charValue(s2[i-1])], col1[j] + delta));
		}
		ans.push_back(col2[m]);
		col1 = col2;
	}

	return ans;
}

pair<string, string> findOptimalCostOptimizedSpace(string s1, string s2, int delta, int alphas[4][4])
{
	int len1 = s1.length();
	int len2 = s2.length();

	if(len1 == 0 && len2 == 0)
	{
		return make_pair("", "");
	}
	if(len1 == 0)
	{
		string p1 = "";
		for(int i=0;i<len2;i++)
		{
			p1 += "_";
		}
		return make_pair(p1, s2);
	}

	if(len2 == 0)
	{
		string p2 = "";
		for(int i=0;i<len1;i++)
		{
			p2 += "_";
		}
		return make_pair(s1, p2);
	}
	if(len1 == 1 && len2 == 1)
	{
		if(s1[0] == s2[0]) return make_pair(s1,s2);
		else
		{
			if(2*delta > alphas[charValue(s1[0])][charValue(s2[0])])
			{
				return make_pair(s1,s2);
			}
			else return make_pair("_"+s1, s2+"_");
		}
	}
	if(len1 == 1 || len2 == 1)
	{
		pair<int, pair<string, string> > ans = findOptimalCost(s1, s2, delta, alphas);
		return ans.second;
	}
	
	int m1 = len1 / 2;

	string leftS1 = s1.substr(0,m1);
	string rightS1 = s1.substr(m1);

	vector<int> s2DividedForLeftS1;
	vector<int> s2DividedForRightS1;
	

	s2DividedForLeftS1 = optimalCostValueOptimisedSpace(leftS1, s2, delta, alphas);


	string s2Reversed = s2;
	string rightS1Reversed = rightS1;
	reverse(s2Reversed.begin(), s2Reversed.end());
	reverse(rightS1Reversed.begin(), rightS1Reversed.end());

	s2DividedForRightS1 = optimalCostValueOptimisedSpace(rightS1Reversed, s2Reversed , delta, alphas);


	reverse(s2DividedForRightS1.begin(), s2DividedForRightS1.end());
	int mini = INT_MAX, ind = -1;
	for(int i=0;i<s2DividedForLeftS1.size();i++)
	{
		if(mini >= s2DividedForLeftS1[i] + s2DividedForRightS1[i])
		{
			mini = s2DividedForLeftS1[i] + s2DividedForRightS1[i];
			ind = i;
		}
	}
	string s2Left = s2.substr(0, ind);
	string s2Right = s2.substr(ind);

	pair<string, string> leftAns = findOptimalCostOptimizedSpace(leftS1, s2Left, delta, alphas);
	pair<string, string> rightAns = findOptimalCostOptimizedSpace(rightS1, s2Right, delta, alphas);

	pair<string, string> finalAns;
	finalAns.first = leftAns.first + rightAns.first;
	finalAns.second = leftAns.second + rightAns.second;

	return finalAns;


}





pair<string, string> preprocessStrings(string s1, string s2, vector<int>& indexesS1, vector<int>& indexesS2)
{

	for(int i=0;i<indexesS1.size();i++)
	{
		string s1_l="", s1_r="";
		int s1Ind=-1, s2Ind=-1;

		if(s1.length() <= indexesS1[i]+1)
			s1_l += s1;
		else{
			s1_l += s1.substr(0,indexesS1[i]+1);
			s1_r +=  s1.substr(indexesS1[i]+1);
		}
		// s1_l += s1.substr(0,indexesS1[i]+1);

		// if(indexesS1[i]+1 > s1.length()-1)
		// 	s1_r +=  s1.substr(indexesS1[i]+1);

		string updatedS1 = "";

		updatedS1 += (s1_l + s1 + s1_r); 

		s1 = updatedS1;
	}

	for(int i=0;i<indexesS2.size();i++)
	{

		string s2_l = "";
		string s2_r = "";
		if(s2.length() <= indexesS2[i]+1)
			s2_l += s2;
		else{
			s2_l = s2.substr(0,indexesS2[i]+1);
			s2_r = s2.substr(indexesS2[i]+1);
		}
		// string s2_r = s2.substr(indexesS2[i]+1);
		string updatedS2 = "";
		updatedS2 += (s2_l + s2 + s2_r);
		s2 = updatedS2;
	}

	return make_pair(s1,s2);

}


void writeOutput(vector<int> iMinCostAdv, string line1, string line2, duration<long long,std::micro> iDuration, float iMemoryConsumed)
{
	ofstream outputFile;
	outputFile.open ("output.txt");
	outputFile << line1 << endl;					//Line1
	outputFile << line2 << endl;				//Line2
	outputFile << iMinCostAdv[iMinCostAdv.size()-1]<<endl;		//Line3 - Cost
	outputFile << setprecision(10) << iDuration.count() / double(1000000) << endl;	//Line4 - Time consumed
	outputFile << setprecision(10) << iMemoryConsumed << endl;						//Line5 - Space consumed - Dummy value (To do - Find actual cost)
	outputFile.close();
	
}
int main(int argc, char* argv[])
{
	auto start = high_resolution_clock::now();

	struct rusage usage;
	int returnValue = getrusage(RUSAGE_SELF, &usage);

 	string fileName = "";
	if(argc == 2)
    {
		fileName = argv[1];	//Input File Name
	}
	if(!fileName.empty())
	{
		ifstream file(fileName);
		string line;
		string s1="", s2 ="";
		vector<int> indexesS1;
		vector<int> indexesS2;
		bool flag = false;
		int num = -1;
		int count = 0;
		while(getline(file, line))
		{
			if(count == 0 && !flag)
			{
				s1 += line;
			}
			else if(line[0] >= 58 && !flag)
			{
				s2 += line;
				flag = true;
			}
			else
			{
				if(!flag){
					if(line=="") continue;
					int val = stoi(line);
					indexesS1.push_back(val);
				}
				else{
					if(num==-1)
						num = 0;
					if(line=="") continue;
					indexesS2.push_back(stoi(line));
					num++;
				}
			}
			count++;
		}


		if(s1.size()>1 && isspace(s1.back()))
		s1 = s1.substr(0, s1.size()-1);

		if(s2.size()>1 && isspace(s2.back()))
		s2 = s2.substr(0, s2.size()-1);
		pair<string, string> alginmentStrings = preprocessStrings(s1, s2, indexesS1, indexesS2);

		string sequence1 = alginmentStrings.first;
		string sequence2 = alginmentStrings.second;
		// cout<<sequence1<<" "<<sequence2<<endl;
		// sequence1 = "TACCCTCAGGGCCGTGGATTTCGTGGTGTTGGACGAACATTTGTGAGCGATGCAGCATCGAAAACTCGTACAGGCCACGATGGTGTTCGCTGTTGGGGGGAAATTTCAGTTTGACACAAGGCGAATAATTGTTTGACTAAACAAAATCATAACGTGCTCCGTGCCGAAAGACGCGAATGCTGTACTTACCATGGGCTGTGCGCATCGCCTCGCGTCCTCCTCGATATGCGCGCATGTAAGGAGTCGGATTATAGTCCATTCAAAGGCACTGTATATGAAGTGCTGTCCGGGAAAGACGCGAATGCTGTACTTACC";
		// sequence2 = "GATCTCGCGACAAACGTCGACTGGAGGTAGGCCGGGAATGGTTGAGACCTATGCTGCCTTTTTGGGTATTTAGTGGCTCGTAGAAAGGGTACTAGACGGTCTGGCCGCATCACAACGACTTATGTTGGAGGTAACGCATTCTTAACCCCACGGAGGAGAGTAACGTGTGCTATTTAGCCGCGTGATCGGTAACCAACTTTCGAGTTCATTATATTAATATACCATTCAAAGGCACTGTATATGAAGTGCTGTCCGGGAAAGACGCGAATGCTGTACTTACCATGGGCTGTGCGCATCGCCTCGCGTCCTCCTCGATAT";

		cout<<sequence1.size()<<" "<<sequence2.size()<<endl;

		int delta = 30;

		int alphas[4][4] = { 	 
			{0,110,48,94},
			{110, 0, 118, 48}, 
			{48,118,0,110}, 
			{94,48,110,0} 
		};

		pair<string, string> minCostDnCnDP = findOptimalCostOptimizedSpace(sequence1, sequence2, delta, alphas);
		vector<int> minCostAdv = optimalCostValueOptimisedSpace(sequence1, sequence2,delta,alphas);
		cout<<minCostAdv[minCostAdv.size()-1]<<endl;
		string firstFiftyX="", firstFiftyY = "";
		string lastFiftyX="", lastFiftyY = "";


		string xAligned = minCostDnCnDP.first;
		string yAligned = minCostDnCnDP.second;
		int lenX = xAligned.length();
		int lenY = yAligned.length();
		if(lenX > 50){
			firstFiftyX = xAligned.substr(0,50);
			lastFiftyX = xAligned.substr(lenX - 50); 
		}
		else
		{
			firstFiftyX = xAligned;
			lastFiftyX = xAligned;
		}

		if(lenY > 50)
		{
			firstFiftyY = yAligned.substr(0,50);
			lastFiftyY = yAligned.substr(lenY - 50); 
		}
		else
		{
			firstFiftyY = yAligned;
			lastFiftyY = yAligned;
		}

		// cout<<firstFiftyX<<" "<<lastFiftyX<<endl;
		// cout<<firstFiftyY<<" "<<lastFiftyY<<endl;
		
		// cout<<minCostDnCnDP.first<<endl;
		// cout<<minCostDnCnDP.second<<endl;


		auto stop = high_resolution_clock::now();

		auto duration = duration_cast<microseconds>(stop - start);
		// cout << duration.count() / double(1000000) << endl;
		
		struct rusage usage2;
        float returnValue2 = getrusage(RUSAGE_SELF, &usage2);
        // cout << float(usage2.ru_maxrss) << endl;
        float memoryConsumed = usage2.ru_maxrss - usage.ru_maxrss;
        memoryConsumed /= float(1024);


        string line1 = firstFiftyX + " " + lastFiftyX;
        string line2 = firstFiftyY + " " + lastFiftyY;

		writeOutput(minCostAdv, line1, line2, duration, memoryConsumed);
		file.close();
	}
}