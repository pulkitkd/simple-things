#include <iostream>
#include <fstream>

using namespace std;

int main()
{
    ifstream fin;
    fin.open("velocity");
    //fout.open("velocity.txt");
    cout<<fin;

    fin.close();
    //fout.close();
    return 0;
}
