import argparse
from gmpy2 import mpz, powmod, jacobi, isqrt
from math import sqrt
from time import time as tm
from tqdm import tqdm as sp

print("riesel.py version 1.0")
res = ""

# Function to compute the Lucas sequence starting value
def find_starting_value(k, n, N):
    if k % 6 == 1 or k % 6 == 5:
        # Calculate u0 using Lucas sequence for N ≡ 7 (mod 24)
        sqrt_3 = sqrt(3)
        term1 = pow(2 + sqrt_3, k)
        term2 = pow(2 - sqrt_3, k)
        u0 = term1 + term2
        return mpz(u0 % N)
    elif k == 3 and (n % 4 == 0 or n % 4 == 3):
        # Special case when k = 3 and n is 0 or 3 mod 4
        return mpz(5778)
    else:
        # General case using Jacobi symbol to find P
        for P in [5, 8, 9, 11]:
            if jacobi(P - 2, N) == 1 and jacobi(P + 2, N) == -1:
                return mpz(powmod(P, k, N))
        return None

# Main function to test for Riesel primes
def riesel(k, n, vfalse=False, vprinti=False, vnores=True):
    global res
    N = mpz(k * (1 << n) - 1)  # Compute N = k * 2^n - 1

    if k >= (1 << n):
        if vprinti:
            res += f"k={k} is greater than 2^n={2**n}, skipping.\n"
        return False
    
    if vprinti:
        res += f"Riesel test for {k}*2^{n}-1\n"
    
    # Find the starting value u0
    u0 = find_starting_value(k, n, N)
    if u0 is None:
        if vfalse:
            res += f"Failed to find u0 for {k}*2^{n}-1!\n"
        return False

    # Initialize u_i = u0 and start the recursion
    u_i = u0
    for i in range(1, n - 1):
        u_i = (u_i * u_i - 2) % N  # Recursively compute u_i
    
    # Check if N divides u_{n-2}
    if u_i == 0:
        res += f"N ({k}*2^{n}-1) is prime with {len(str(N))} digits.\n"
        return True
    else:
        if vfalse:
            res += f"N ({k}*2^{n}-1) is not prime.\n"
        return False

# Main program loop to parse arguments and run the test
def main():
    global res
    parser = argparse.ArgumentParser(description="Riesel Prime Test")
    parser.add_argument("--min-k", type=int, default=1, help="Minimum k value (default: 1)")
    parser.add_argument("--max-k", type=int, default=10, help="Maximum k value (default: 10)")
    parser.add_argument("n", type=int, help="The n value in k*2^n-1")
    parser.add_argument("--max-a", type=int, default=10, help="Maximum value of a to test (default: 10)")
    parser.add_argument("--vfalse", action="store_true", help="Verbose output when test fails")
    parser.add_argument("--vprinti", action="store_true", help="Verbose output including initial test information")
    parser.add_argument("--vnores", action="store_true", help="Verbose output when no result is conclusive")
    
    args = parser.parse_args()
    have_prime = False
    t = tm()

    # Iterate over odd k values from min-k to max-k
    for k in sp(range(args.min_k // 2 * 2 + 1, args.max_k + 1, 2)):
        try:
            result = riesel(k, args.n, args.vfalse, args.vprinti, args.vnores)
        except Exception as e:
            print(f"Error during testing: {e}")
            print(res)
            return
        
        if result:
            have_prime = True

    if not have_prime:
        res += "No prime found.\n"
    
    res += f"Time: {tm() - t:.6f} seconds\n"
    print(res)

if __name__ == "__main__":
    main()
