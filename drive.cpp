#include "PML.h"
int main(){
	//PML pml(5, 5, 100, 1, "C:\\Users\\Leo\\Desktop\\1.txt", "C:\\Users\\Leo\\Desktop\\2.txt");
	PML pml(100, 100, 1000, 1, "1.txt", "2.txt");
	pml.Solve();
	//system("pause");
	return 1;
}