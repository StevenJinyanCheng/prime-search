#include <iostream>
#include <sstream>
#include <string>
#include "C:\ProgramData\anaconda3\Library\include\gmp.h"
#include <ctime>
#include <vector>
#include <stdexcept>

// Global result string
std::string res;

// Function to perform Proth prime test
bool proth(int k, int n, int max_a = 10, bool vfalse = false, bool vprinti = false, bool vnores = true) {
    bool verbose = true;
    if (n < 32 && k > (1ll << n)) {
        if (verbose) {
            res += "k=" + std::to_string(k) + " is greater than 2^n=" + std::to_string(1 << n) + ", skipping.\n";
        }
        return false;
    }
    if (verbose && vprinti) {
        res += "Proth test " + std::to_string(k) + "*2^" + std::to_string(n) + "+1\n";
    }
    mpz_t p_12, p_1, p, powres;
    mpz_init(p_12);
    mpz_init(p_1);
    mpz_init(p);
    mpz_init(powres);

    mpz_set_ui(p_12, k);
    mpz_mul_2exp(p_12, p_12, n - 1);

    mpz_set_ui(p_1, k);
    mpz_mul_2exp(p_1, p_1, n);
    mpz_add_ui(p, p_1, 1);

    for (int i = 2; i <= max_a; i++) {
        mpz_t tmp;
        mpz_init_set_ui(tmp, i);
        mpz_powm(powres, tmp, p_12, p);  // powres = tmp^p_12 mod p
        mpz_clear(tmp);
        if (mpz_cmp(powres, p_1) == 0) {
            if (verbose) {
                res += "Find " + std::to_string(i) + "^((p-1)/2)=-1 mod p, p (" + std::to_string(k) + "*2^" + std::to_string(n) + "+1) is a prime with " + std::to_string(mpz_sizeinbase(p, 10)) + " digits.\n";
            }
            mpz_clears(p_12, p_1, p, powres, NULL);
            return true;
        }
        if (mpz_cmp_ui(powres, 1) != 0) {
            if (verbose && vfalse) {
                res += "p (" + std::to_string(k) + "*2^" + std::to_string(n) + "+1) is not a prime because find " + std::to_string(i) + "^((p-1)/2) != 1 and != -1 mod p\n";
            }
            mpz_clears(p_12, p_1, p, powres, NULL);
            return false;
        }
    }

    if (verbose && vnores) {
        res += "Proth test inconclusive for " + std::to_string(k) + "*2^" + std::to_string(n) + "+1\n";
    }

    mpz_clears(p_12, p_1, p, powres, NULL);
    return false;
}

int main(int argc, char* argv[]) {
    int min_k = 1, max_k = 10, n = 0, max_a = 10;
    bool vfalse = false, vprinti = false, vnores = false;

    // Parsing command line arguments
    for (int i = 1; i < argc; i++) {
        std::string arg(argv[i]);
        if (arg == "--min-k" && i + 1 < argc) {
            min_k = std::stoi(argv[++i]);
        } else if (arg == "--max-k" && i + 1 < argc) {
            max_k = std::stoi(argv[++i]);
        } else if (arg == "--max-a" && i + 1 < argc) {
            max_a = std::stoi(argv[++i]);
        } else if (arg == "--vfalse") {
            vfalse = true;
        } else if (arg == "--vprinti") {
            vprinti = true;
        } else if (arg == "--vnores") {
            vnores = true;
        } else if (n == 0) {  // Ensure n is set only once and correctly
            n = std::stoi(argv[i]);
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            return 1;
        }
    }

    // Debugging output
    std::cout << "Parsed arguments:\n";
    std::cout << "min_k: " << min_k << "\n";
    std::cout << "max_k: " << max_k << "\n";
    std::cout << "n: " << n << "\n";
    std::cout << "max_a: " << max_a << "\n";
    std::cout << "vfalse: " << vfalse << "\n";
    std::cout << "vprinti: " << vprinti << "\n";
    std::cout << "vnores: " << vnores << "\n";

    if (n <= 0) {
        std::cerr << "n must be greater than 0.\n";
        return 1;
    }

    bool havep = false;
    clock_t start = clock();

    try {
        int cnt = 1, total = 0;
        for (int k = (min_k / 2) * 2 + 1; k <= max_k; k += 2, total++);
        for (int k = (min_k / 2) * 2 + 1; k <= max_k; k += 2, cnt++) {  // Ensure k is odd
            if (proth(k, n, max_a, vfalse, vprinti, vnores)) {
                havep = true;
            }
            printf("%d/%d %.2f%%\r", cnt, total, (cnt*100.0f/total));
        }
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << "\n";
        std::cerr << res;
        return 1;
    }

    if (!havep) {
        res += "There is no prime found.\n";
    }

    double elapsed = static_cast<double>(clock() - start) / CLOCKS_PER_SEC;
    res += "Time: " + std::to_string(elapsed) + " seconds\n";
    std::cout << std::endl << res;

    return 0;
}
