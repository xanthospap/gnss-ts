#include <iostream>
#include "flags.hpp"

int main ()
{
    // initialize a flag; all fields should be '0'
    flag f;
    std::cout << "\nDefault initialized flag: ";
    f.debug_print();
    std::cout << "\n";

    // set the outlier flag to  'true'
    f.set_outlier();
    // ask the flag if it is an outlier
    if ( f.is_outlier() ) {
      std::cout << "Outlier bit set!";
    } else {
      std::cout << "Not an outlier!";
    }
    std::cout << "\n";
    // print the flags
    std::cout << "\nOutlier flag set        : ";
    f.debug_print();
    std::cout << "\n";

    // unset the outlier flag
    f.unset_outlier();
    // ask the flag if it is an outlier
    if ( f.is_outlier() ) {
      std::cout << "Outlier bit set!";
    } else {
      std::cout << "Not an outlier!";
    }
    std::cout << "\n";
    // print the flags
    std::cout << "\nOutlier flag unset      : ";
    f.debug_print();

    // lets set two flags, both an earthquake and a jump
    f.set_outlier();
    std::cout << "\nFlag set                : ";
    f.debug_print();
    std::cout << "\n";
    f.set_skip();
    std::cout << "\nFlag set                : ";
    f.debug_print();
    std::cout << "\n";
    f.set_erthq();
    std::cout << "\nFlag set                : ";
    f.debug_print();
    std::cout << "\n";
    f.set_velchg();
    std::cout << "\nFlag set                : ";
    f.debug_print();
    std::cout << "\n";
    f.set_jump();
    std::cout << "\nFlag set                : ";
    f.debug_print();
    std::cout << "\n";
    
    // ask the flag if it is an outlier
    if ( f.is_outlier() ) {
      std::cout << "Outlier bit set!";
    } else {
      std::cout << "Not an outlier!";
    }
    std::cout << "\n";
    // ask if it is a skip
    if ( f.is_skip() ) {
      std::cout << "Skip bit set!";
    } else {
      std::cout << "Not a skip!";
    }
    std::cout << "\n";
    // ask the flag if it is an earthquake
    if ( f.is_erthq() ) {
      std::cout << "Earthquake bit set!";
    } else {
      std::cout << "Not an earthquake!";
    }
    std::cout << "\n";
    // ask if it is a velocity change
    if ( f.is_velchg() ) {
      std::cout << "Vel. change bit set!";
    } else {
      std::cout << "Not a velocity change!";
    }
    std::cout << "\n";
    // ask the flag if it is a jump
    if ( f.is_jump() ) {
      std::cout << "Jump bit set!";
    } else {
      std::cout << "Not a jump!";
    }

    std::cout << "\n";
    return 0;
}
