#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <cstring>
#include <cctype>
#include <cstdlib>
#include <set>
#include <vector>
#include <map>

using namespace std;

int globe_coverThresh;
double globe_mutThresh, globe_insThresh, globe_delThresh; 
double globe_samp_mutThresh;

struct pileup
{
 string chr;
 int pos;
 char base;
 int cover;
 string readStr;
 string quality;
 void clear();
};

class midPerc
{
 public: 
 int coverMut;
 map<char, double> bases;

 int coverDel;
 double delPerc;

 int  coverIns;
 int numOfIns; //2
 double insPerc;
 double dorminInsPerc;
 string insStr; //"AC"
 
 public:
 midPerc()
 {
  bases['A'] = 0;
  bases['C'] = 0;
  bases['G'] = 0;
  bases['T'] = 0;
 };
 bool isMut(vector<string> & resStr, const midPerc & b) const;
 bool isDel(vector<string> & resStr, const midPerc & b) const;
 bool isIns(vector<string> & resStr, const midPerc & b) const;
 bool detectMID(vector<string> & resStr, const midPerc & b) const;
 void clear();
};

//parse +3ACG    -2TT strings, obtain ACG & TT
int parseInDelStr(string mapStr, int i, string & ins);
//use the content of pileup to calculate mut, ins, del percentages
void pileup2midPerc(midPerc & mid, const pileup & pile);
//parse a line in a pileup file into pile object.
//if the line has < 6 fields, return false
bool line2pile(string line, pileup & pile);
string getOpFileName(char *, char *);
bool setChromOrder(char *, map<string, int> &);
void detectSNP(ifstream & input1, ifstream & input2, ofstream & output, const map<string, int> & );
int compareTwoLines(string line1, string line2, pileup & pi1, pileup & pi2, const  map<string, int> & chrOrder);

void pileup::clear()
{
 chr = "";
 pos = 0;
 base = '0';
 cover = 0;
 readStr = "";
 quality = "";
}

void midPerc::clear()
{
 coverMut = 0;
 bases['A'] = 0;
 bases['C'] = 0;
 bases['G'] = 0;
 bases['T'] = 0;


 coverDel = 0;
 delPerc = 0;

 coverIns = 0;
 numOfIns = 0; //2
 insPerc = 0;
 dorminInsPerc = 0;
 insStr = ""; //"AC"
 
}

bool midPerc::isMut(vector<string> & resStr, const midPerc & b) const 
{
 resStr.clear();
 if(coverMut < globe_coverThresh || b.coverMut < globe_coverThresh)
   return false;
 
 map<char, double>::const_iterator iter;
 double t1 = 0, t2 = 0, control_c2_perc, samp_c2_perc;
 char c1 = '#', c2 = '#';
 string str, controlGeno = "", sampGeno = "";
// ****** coverMut == 6
// c1 or c2 undefined
 for(iter = bases.begin(); iter != bases.end(); iter++)
 {
  if(iter->second > t1)
  {
   t1 = iter->second;
   c1 = iter->first;
  }
 }
 for(iter = bases.begin(); iter != bases.end(); iter++)
 {
  if( iter->second > 0.3)
    controlGeno += iter->first;  //aacccc -> AC; aa*** -> A; a*** -> false; **** -> false;
 }
 for(iter = b.bases.begin(); iter != b.bases.end(); iter++)
 {
  if(iter->second > t2)
  {
   t2 = iter->second;
   c2 = iter->first;
  }
  if( iter->second > 0.3)
     sampGeno += iter->first;
 }
 if(c1 == '#' || c2 == '#' || controlGeno == "")
   return false;
// 1. sample's major allele c2 >= 0.90
// 2. in control, c2's perc <= 0.05
// therefore, in control, the major allele c1 must be >= 0.95/3 && c1 != c2
// c1 -> c2 in the output 
 for(iter = bases.begin(); iter != bases.end(); iter++)
 {
  if(iter->first == c2)
  {
   control_c2_perc = iter->second;
   break;
  }
 }
 if(t2 >= globe_samp_mutThresh &&  control_c2_perc  <= globe_mutThresh)
 {
  str = controlGeno;
  str = str + "->" + c2;
  resStr.push_back("M");
  resStr.push_back(str);
  return true;
 }
 for(iter = b.bases.begin(); iter != b.bases.end(); iter++)
 {
  if(iter->first == c1)
  {
   samp_c2_perc = iter->second;
   break;
  }
 }
 if(t1 >= globe_samp_mutThresh &&  samp_c2_perc  <= globe_mutThresh)
 {
  str = c1;
  str = str + "->" + sampGeno;
  resStr.push_back("M");
  resStr.push_back(str);
  return true;
 }
 return false;
}

bool midPerc::isDel(vector<string> & resStr, const midPerc & b) const
{
 resStr.clear();
 if(coverDel < globe_coverThresh || b.coverDel < globe_coverThresh)
   return false;
 string str;
 if(delPerc >= globe_delThresh && b.delPerc <= 0.1)
 {
  resStr.push_back("I");
  resStr.push_back("1D");
  return true;
 }
 if(delPerc <= 0.1 && b.delPerc  >= globe_delThresh)
 {
  resStr.push_back("D");
  resStr.push_back("2D");
  return true;
 }

 return false;
}

bool midPerc::isIns(vector<string> & resStr, const midPerc & b) const
{
 resStr.clear();
 if(coverIns < globe_coverThresh || b.coverIns < globe_coverThresh)
   return false;
 string str, s1, s2;
 stringstream ss;

 if(insPerc >= globe_insThresh || b.insPerc >= globe_insThresh)
 {
  if(b.insPerc <= 0.1)
  {
   ss.str("");
   ss.clear();
   ss << numOfIns;
   str = "+";
   str = str + ss.str() + insStr; //"+3AGC"
   str = str + "->+0";
   resStr.push_back("D");
   resStr.push_back(str);
   return true;
  }
  if(insPerc <= 0.1) // b have insertion
  {
   ss.str("");
   ss.clear();
   ss << b.numOfIns;
   str = "+0->";
   str = str + "+" +  ss.str() + b.insStr; //"+0->+3AGC"
   resStr.push_back("I");
   resStr.push_back(str);
   return true;
  }
  if(insPerc >= globe_insThresh && b.insPerc >= globe_insThresh && dorminInsPerc  >= globe_insThresh && b.dorminInsPerc  >= globe_insThresh)
  {
    s1 = "+";
    ss.str("");
    ss.clear();
    ss << numOfIns;
    s1 = s1 + ss.str() + insStr; //"+3AGC"
    s2 = "+";
    ss.str("");
    ss.clear();
    ss << b.numOfIns;
    s2 = s2 + ss.str() + b.insStr; //"+3AGC"
    
    if(numOfIns == b.numOfIns) 
    {
      if(insStr != b.insStr)
      {
       resStr.push_back("M");
       resStr.push_back(s1 + "->" + s2);
       return true;
      }
    }
    if(numOfIns > b.numOfIns)
    {
     resStr.push_back("D");
     resStr.push_back(s1 + "->" + s2);
     return true;
    }
    if(numOfIns < b.numOfIns)
    {
     resStr.push_back("I");
     resStr.push_back(s1 + "->" + s2);
     return true;
    }
  }
 }
 return false;
}

bool midPerc::detectMID(vector<string> & resStr, const midPerc & b) const
{
 bool m, i, d;
 string str1 = "", str2 = "";
 vector<string> s1, s2, s3;
 m = isMut(s1, b);
 i = isIns(s2, b);
 d = isDel(s3, b);
 if(m || i || d) //at least one is true
 {
  if(m)
  {
   str1 = s1[0]; //"M"
   str2 = s1[1] + ";"; // A->C
  }
  if(i)
  {
   str1 += s2[0]; // "MI"
   str2 += s2[1]  + ";"; //A->C;+1A->+1C
  }
  if(d)
  {
   str1 += s3[0]; // "MID"
   str2 += s3[1] + ";"; //A->C;+1A->+1C;2D
  }
  resStr.clear();
  resStr.push_back(str1);
  resStr.push_back(str2);
  return true;
 }
 return false;
}

int parseInDelStr(string mapStr, int i, string & ins)  //+20TTTTTTTTTTTT....TG; -3AAC
{
   int begin ;
   i++; //the first digit
   string strIn(1, mapStr[i]);
   i++;
   while(isdigit(mapStr[i]))
   {
    strIn += mapStr[i];
    i++;
   }
   //exit while loop, i points to the first 'acgtn'
   begin = i;
   int numIn = atoi(strIn.c_str());
   i = i - 1 + numIn; //i points to the last 'agctn'
   int length = i + 1 - begin;
   ins = mapStr.substr(begin, length);
   return i;
}

void pileup2midPerc(midPerc & mid, const pileup & pile)
{
 char letter, refBase;
 string mapStr = "", insStr = "", acg = "ACGT"; 
 vector<string> insList;
 int numIn = 0, numDel = 0;
 int numMatches = 0, mutCover = pile.cover, insCover = pile.cover, delCover = pile.cover;
 map<char, int> baseCount;
 map<string, int> insDormMap;
 bool res = false;

 for(int i = 0; i < acg.length(); i++)
   baseCount[acg[i]] = 0;
 refBase = toupper(pile.base);
 for(int i = 0; i < pile.readStr.length(); i++)
   mapStr += toupper(pile.readStr[i]);

 for(int i = 0; i < mapStr.length(); i++)
 {
  letter = mapStr[i];
  if(letter == '^')
  {
   i++; //skip the char following ^
   delCover--;
  }
  else if (letter == '$')
  {
   insCover--;
   delCover--;
  }
  else if (letter == '+' ) //+20TTTTTTTTTTTT....TG
  {
   i = parseInDelStr(mapStr, i, insStr); //i points to the last 'agctn'
   insList.push_back(insStr);
   numIn++;
  }
  else if (letter == '-' )
  {
   i = parseInDelStr(mapStr, i, insStr); //i points to the last 'agctn'
  }
  else if (letter == '*' )
  {
   numDel++;
  }
  else if (letter == ',' || letter == '.')
  {
   numMatches++;
  }
  else if (letter == 'N' || letter == 'n')
  {
   mutCover--;
  }
  else if(acg.find(letter) != string::npos) //ACGT acgt
  {
   baseCount[letter]++;
  }
  else; //other symbol
 }
 mid.coverMut = mutCover;
 if(mid.coverMut >= globe_coverThresh)
 {
  //populate bases
  baseCount[refBase] = numMatches;
  for(int i = 0; i < acg.length(); i++)
    mid.bases[acg[i]] = baseCount[acg[i]] / (double) mid.coverMut; 
 }
 mid.coverDel = delCover;
 if(mid.coverDel  >= globe_coverThresh)
 {
  mid.delPerc = numDel / (double) mid.coverDel;
 }
 mid.coverIns = insCover;
 if(mid.coverIns >= globe_coverThresh)
 {
  mid.insPerc = numIn / (double)  mid.coverIns;
  if(numIn > 0)
  {
   for(int i = 0; i < insList.size(); i++)
   {
    if(insDormMap.count(insList[i]) == 0)
    {
     insDormMap[insList[i] ] = 1;
    }
    else
      insDormMap[insList[i] ]++;
   }
   int tempInt = 0;
   string tempStr;
   map<string, int>::const_iterator iter;
   for( iter = insDormMap.begin(); iter != insDormMap.end(); iter++)
   {
     if(iter->second > tempInt)
     {
       tempInt = iter->second;
       tempStr = iter->first;
     }
   }
   mid.insStr = tempStr;
   mid.numOfIns = tempStr.length();
   mid.dorminInsPerc = tempInt / (double) numIn; 
  }
 }
}

bool line2pile(string line, pileup & pile)
{
 int numTab = 0;
 vector<int> tabPos;

 for(int i = 0; i < line.length(); i++)
 {
  if(line[i] == '\t')
  {
   tabPos.push_back(i);
  }
 }
 if(tabPos.size() < 5)
   return false;
 int len, coverage;
 char refBase;
 string mapStr, str;
 
 //chrosome 1st field
 len = tabPos[0] ;
 pile.chr = line.substr(0, len);
 //chr position 2nd field, between 1st and 2nd tab
 len = tabPos[1] - tabPos[0] - 1;
 str = line.substr(tabPos[0] + 1, len);
 pile.pos = atoi(str.c_str());
 //ref base is following the 2nd tab
 refBase = line[tabPos[1] + 1 ];
 pile.base = refBase;
 //coverage between 3rd tab and 4th tab
 len = tabPos[3] - tabPos[2] - 1;
 str = line.substr(tabPos[2] + 1, len);
 coverage = atoi(str.c_str());
 pile.cover = coverage;
 //length should be 4th tab and 5th tab
 len = tabPos[4] - tabPos[3] - 1;
 mapStr = line.substr(tabPos[3] + 1, len);
 pile.readStr = mapStr;
 return true;
}

string getOpFileName(char *f1, char *f2)
{
 string s1 = f1, s2 = f2;
 string res;
 int i1 = s1.find_last_of('/');
 int i2 = s2.find_last_of('/');
 if(i1 != string::npos)
   s1 = s1.substr(i1 + 1);
 if(i2 != string::npos)
   s2 = s2.substr(i2 + 1);
 res = s1 + "_" + s2 + ".snp";
 return res;
}

bool setChromOrder(char *a, map<string, int> & b)
{
 ifstream input;
 string chr;
 int order = 0; 

 input.open(a, ios::in);
 if(!input)
   return false;
 getline(input, chr);
 while(!input.eof())
 {
  if(chr.length() > 0)
  {
   b[chr] = order;
   order++;
  }
  getline(input, chr);
 }
 return true;
}

void detectSNP(ifstream & input1, ifstream & input2, ofstream & output, const map<string, int> & chrOrder)
{
 string line1, line2;
 pileup pi1, pi2;
 midPerc m1, m2;
 int compRes;
 vector<string> midCode;

 getline(input1, line1);
 getline(input2, line2);
 while(!input1.eof() && !input2.eof())
 {
  pi1.clear();
  pi2.clear();
  m1.clear();
  m2.clear();
  midCode.clear();
  compRes = compareTwoLines(line1, line2, pi1, pi2, chrOrder);
  if(compRes < 0)
    getline(input1, line1);
  else if(compRes > 0)
    getline(input2, line2);
  else
    {
     pileup2midPerc(m1, pi1);
     pileup2midPerc(m2, pi2);
     if(m1.detectMID(midCode , m2) == true)
     {
      char t = '\t';
      output << midCode[0] << t << midCode[1] << t << pi1.chr << t << pi1.pos << t <<  pi1.base;
      output << t << pi1.cover << t << pi1.readStr;
      output << t << pi2.cover << t << pi2.readStr << endl;
     }
     getline(input1, line1);
     getline(input2, line2);
    }
 }
}

int compareTwoLines(string line1, string line2, pileup & pi1, pileup & pi2, const  map<string, int> & chrOrder)
{
 map<string, int>::const_iterator it1, it2;
 if(line2pile(line1, pi1) == false)
   return -1; //e.g. an empty line in pileup1 file, read next line
 // pi1 has 6 fields
 if(line2pile(line2, pi2) == false)
   return 1;
 //pi2 has 6 fields
 if(pi1.chr == pi2.chr)
 {
  return (pi1.pos - pi2.pos);  
 }
 //pi1, pi2 in different chromoses
 //pi1.chr order < pi2.chr order , 1 read next line, return < 0
 //pi1.chr order > pi2.chr order, 2 read next line, return > 0
 //return (chrOrder[pi1.chr] - chrOrder[pi2.chr]);
 it1 = chrOrder.find(pi1.chr);
 it2 = chrOrder.find(pi2.chr);
 if(it1 == chrOrder.end()) //file 1 line's chromosome is not in the chr-order-file, read next line
   return -1;
 if(it2 == chrOrder.end())
   return 1;
 //both chr in the chr-order-file
 return it1->second - it2->second ;
}

int main(int argc, char *argv[])
{
 if(argc != 8)
 {
  cout << argv[0] << "  pileup-file-1  pileup-file-2  same-allele-max-perc(0.05) 	gapCutoff(0.8) chromosome-order-file  coverage_cutoff(5) sampMutCutoff(0.85) \n";
  cerr << "sample must be homozyg, e.g. 85%, 90%" << endl;
  cerr << "control's allele, which is same as sample's dorminant allele, must be less than, say 5%" << endl;
  cerr << "control's major allele must >= 0.95/3, and it must be different than sample's major allele" << endl;
  return 1;
 }

 string opFile;
 map<string, int> chromOrder;
 ifstream input1, input2;
 ofstream output;
 
 globe_mutThresh = atof(argv[3]);
 globe_insThresh = globe_delThresh = atof(argv[4]);
 globe_coverThresh = atoi(argv[6]);
 globe_samp_mutThresh = atof(argv[7]);
 input1.open(argv[1], ios::in);
 if(!input1)
   {
    cerr << argv[1] << " pileup file cannot be opened for reading!" << endl;
    return 1;
   }

 input2.open(argv[2], ios::in);
 if(!input2)
   {
    cerr << argv[2] << " pileup file cannot be opened for reading!" << endl;
    return 1;
   }
 
 if(!setChromOrder(argv[5], chromOrder))
 {
  cerr << argv[5] << ": error in chromosome order file." << endl;
  return 1;
 }

 opFile = getOpFileName(argv[1], argv[2]);
 output.open(opFile.c_str(), ios::out);
 if(!output)
   {
    cerr << opFile  << " cannot be opened for writing" << endl;
    return 1;
   }
 
//if two lines match to same chr-coord,
// both files read the next line
//no match
// if two lines point to same chr, line2 coord < line1 coord
//only read next line in file2
//
 detectSNP(input1, input2, output, chromOrder);

 input1.close();
 input2.close();
 output.close();
 return 0;
}
