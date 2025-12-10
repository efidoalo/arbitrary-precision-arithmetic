#include "precision.h"
#include <stdlib.h>
#include <stdio.h>

int main()
{
	struct num *n1 = new_num(3,2);
	n1->binary_digits[0] = 0;
	n1->binary_digits[1] = 1;
	n1->binary_digits[2] = 143;
	struct num *n2 = new_num(1,0);
	n2->binary_digits[0] = 2;
	printf("n2: %s  n1: %s\n", num_to_char(n2), num_to_char(n1));
	struct num *n3 = div_num(n2, n1);
	printf("n2 / n1: %s\n", num_to_char(n3));
	unsigned long int n = 1234;
	struct num *n4 = 0;
	set_num_int(&n4, n, 0);
	printf("\n%s\n", num_to_char(n4));
}
