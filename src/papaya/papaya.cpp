#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;

class Contig {
public:
  Contig(string p_name, int p_size);
  void operator+=(int toAdd);
  bool operator<(Contig &A);
  bool operator>(Contig &A);
  int getSize();
  string getName();
private:
  string name;
  int size;
};

Contig::Contig(string p_name, int p_size) : name(p_name), size(p_size) { }

void Contig::operator+=(int toAdd) {
  size += toAdd;
}

bool Contig::operator<(Contig &A) {
  return (size < A.getSize());
}

bool Contig::operator>(Contig &A) {
  return (size > A.getSize());
}

int Contig::getSize() {
  return size;
}

string Contig::getName() {
  return name;
}

int main(int argc, char ** argv){
	if(argc < 2 || argc > 3){
		cout << "papaya [Fasta file] (int) > output" << endl;
		exit(0);
	}
	string input(argv[1]);
	string cline;
	ifstream in(input);
	if(!in){
		cout << "Problem opening .fa(sta) file" << endl;
		return 0;
	}
  int n(0);
  if (argc > 2) {
    string nin(argv[2]);
    n = strtol(nin.c_str(), NULL, 10);
  }

	vector<Contig> ctgs;
  int cn(0);

	while(not in.eof()){
		getline(in,cline);
    if (cline[0] == '>') {
      ctgs.push_back(Contig(cline.substr(1, cline.find(" ")), 0));
      cn++;
    } else {
      ctgs[cn-1] += cline.size();
    }
	}

  sort(ctgs.begin(),ctgs.end());
  if (n == 0) {
    n = ctgs.size();
  }
  int p(0);
  for (int i(ctgs.size()-1); i != 0; i--) {
    if (p < n) {
      cout << ctgs[i].getName() << "\t0\t" << ctgs[i].getSize() << endl;
      p++;
    } else {
      break;
    }
  }
}
