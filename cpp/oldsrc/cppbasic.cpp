#include "cppbasic.hpp"
#include <algorithm>

typedef std::vector<std::string> stringVec;
using   std::vector;
using   std::string;

stringVec lineToStringVector(const std::string& s, const char c)
{
	stringVec ret;
	if (!s.size())
		return ret;
	string newString;
	string scopy (c+s);
	unsigned int i,size(scopy.size()-1);
	for (i=0;i<size;i++) {
		if (scopy[i]==c&&scopy[i+1]!=c) {
			newString=scopy[i+1];
			i+=2;
			while (scopy[i]&&scopy[i]!=c)
				newString+=scopy[i++];
			ret.push_back(newString);
			i--;
		}
	}
	return ret;
}

string strip(const string& str, const char c)
{
	string s (str);
	while (*s.begin()==c)
		s.erase(s.begin());
	while ( *(s.end()-1)==c )
		s.erase(s.end()-1);
	return s;
}

std::string for2cppfloat(const std::string& s)
{
  std::string str (s);
  replace(str.begin(),str.end(),'D','e');
  replace(str.begin(),str.end(),'d','e');
  return str;
}
