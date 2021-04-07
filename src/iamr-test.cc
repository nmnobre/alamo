#include "Util/Util.H"

#include "NavierStokes.H"

// 
// THIS FILE SHOULD EVENTUALLY BE DELETED
//
// The purpose of this file is to test including and linking against IAMR libraries.
// 

int main()
{
    Util::Initialize();

    NavierStokes myns;

    Util::Finalize();
}