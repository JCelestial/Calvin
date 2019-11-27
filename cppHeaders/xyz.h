/*
Header file for any input files of the *.xyz extension

This object will improve upon the previous version
of CALVIN where it uses a sepcific algorithm to
detect the termination sequence of the array but
instead, using a different type of data structure:
The Linked List. This will eliminate the need for the
Shell Script

This will give the program more flexibility in terms
of scalability of code.

*/

#ifndef xyz_h
#define xyz_h

class xyz {
private:

//########## ATTRIBUTES ###########

	typedef struct xyznode {
		// All LL properties
		char atomid[3];
		double xcoor;
		double ycoor;
		double zcoor;
		// Pointer to next element in pdb
		xyznode* next;
	}* xyzPtr;

	xyzPtr head;
	xyzPtr current;
	xyzPtr temporary;

//########## METHODS ###########

public:
//functions

	xyz(); //Linked List constructor
	void AddNodeXyz(pdbPtr *addNode);
	void DeleteNodeXyz(pdbPtr *delNode);
	void PrintTraj(); // prints out elements


};


#endif