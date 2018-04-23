#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

// Armadillo documentation is available at:
// http://arma.sourceforge.net/docs.html

int
main(int argc, char** argv)
{
    mat A(2,3);

    A.set_size(4,5);
    A << 0.165300 << 0.454037 << 0.995795 << 0.124098 << 0.047084 << endr
        << 0.688782 << 0.036549 << 0.552848 << 0.937664 << 0.866401 << endr
        << 0.348740 << 0.479388 << 0.506228 << 0.145673 << 0.491547 << endr
        << 0.148678 << 0.682258 << 0.571154 << 0.874724 << 0.444632 << endr
        << 0.245726 << 0.595218 << 0.409327 << 0.367827 << 0.385736 << endr;
      
    rowvec a = A.row(0);
    a.print("A(0):");

    a.save("exercise1.txt", raw_ascii);

    vec b(5), c(5), d(6);

    d.print("bleh");

    //dot(b, c) // is a thing (documentation for armadillo is back online )
    // as of April 7, 2018

    // Notes:   vec is automatically a column vector and therefore core dumps
    //          when you try to print out a row of a matrix saved into a vec

    cout << "b.n_rows: " << b.n_rows << endl; // output: 5 

  return 0;
}