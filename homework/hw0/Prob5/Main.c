#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

void approximate_pi(double threshold) {
    double diff, pi_appx = 0, pi_true = acos(-1.);
    int n = 0;
    while(n < 50) {
        diff = fabs(pi_appx - pi_true);
        if (diff > threshold) {
            pi_appx = pi_appx +
                pow(16, -n)*(4./(8*n + 1) - 2./(8*n + 4) - 1./(8*n + 5) - 1./(8*n + 6));
        } else {
            break;
        }
        n++;
    }

    printf("Number of iterations: %d\n", n);
    printf("Threshold value: %e\n", threshold);
    printf("PI True:         %1.16f \n", pi_true);
    printf("PI Approximated: %1.16f \n", pi_appx);
    printf("\n\n");
}


int main(int argc, char **argv) {
    double threshold = 1.e-4;
    approximate_pi(threshold);
    threshold = 1.e-8;
    approximate_pi(threshold);
    threshold = 1.e-12;
    approximate_pi(threshold);
    threshold = 1.e-16;
    approximate_pi(threshold);

    return 0;
}
