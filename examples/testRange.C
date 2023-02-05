#include <iostream>
#include <range.h>

// load first range library
//
//   root [0] gSystem->Load("/usr/local/lib/librange")
//
// then compile
//
//   root [1] .L testRange.C++
//
// and run
//
//   root [2] testRange()

void testRange() {

  Double_t err;

  cout << passage(1,2,4,0,14,28,100.,15.,&err) << endl;
  cout << passage(1,2,4,0,14,28,120.,15.,&err) << endl;
  cout << passage(1,2,4,0,14,28,140.,15.,&err) << endl;
  cout << passage(1,3,6,0,14,28,160.,15.,&err) << endl;
  cout << passage(1,3,7,0,14,28,280.,150.,&err) << endl;

}
