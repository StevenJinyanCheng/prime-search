import argparse
from gmpy2 import mpz, powmod, gcd
from time import time as tm
from tqdm import tqdm as sp
def proth_test(k, n, max_a=10):
    res = ""
    if k > (1 << n):
        res += (f"k={k} is greater than 2^n={2**n}, skipping.\n")
        return False
    p_12 = mpz(k << (n - 1))
    p_1 = mpz(p_12 << 1)
    p = mpz(p_1 + 1)
    for i in range(2, max_a + 1):
        powres = powmod(i, p_12, p)
        if powres == p_1:
            res += (f"Find {i}^((p-1)/2)=-1 mod p, p ({k}*2^{n}+1) is a prime with {len(str(p))} digits.\n")
            return True, res
        if powres != 1:
            res += (f"p ({k}*2^{n}+1) is not a prime because find {i}^((p-1)/2) != 1 and != -1 mod p\n")
            return False, res
    res += (f"Proth test inconclusive for {k}*2^{n}+1\n")
    return False, res
def trial_division(p, ub=100, pp=True):
    """Perform trial division up to bound"""
    l = range(2, ub+1)
    if pp:
        l = sp(l)
    for i in l:
        if p % i == 0:
            return False, f"p is divisible by {i}, hence not a prime.\n"
    return True, f"p passed trial division up to {ub}.\n"
def p_minus_1_factorization(p, B=1000):
    """Perform the p-1 factorization method"""
    a = mpz(2)  # Usually a small number like 2 is used for p-1 factorization
    q = mpz(2)
    exp = mpz(1)
    
    # Multiply exp by all primes q up to B, and the highest power of q that divides B
    while q <= B:
        q_pow = q
        while q_pow <= B:
            exp *= q_pow
            q_pow *= q
        q += 1
    
    gcd_val = gcd(powmod(a, exp, p) - 1, p)
    if 1 < gcd_val < p:
        return gcd_val, f"p-1 factorization found a non-trivial factor {gcd_val} of p.\n"
    return None, f"p-1 factorization did not find any factors for p.\n"
def main():
    parser = argparse.ArgumentParser(description="Show details of a k*2**n+1 with k<2**n")
    parser.add_argument("-n", type=int, help="The n value", required=True)
    parser.add_argument("-k", type=int, help="The k value", required=True)
    parser.add_argument("-t", type=int, help="The max trial division value", required=False, default=2000000)
    parser.add_argument("-f", action='store_true', help="factor or not?")
    parser.add_argument("-B", type=int, help="factor bound", required=False, default=10)
    args = parser.parse_args()
    p = 2**args.n*args.k+1
    p = mpz(p)
    res, rstr = trial_division(p, ub=1000, pp=False)
    print(rstr)
    if res:
        res, rstr = trial_division(p, ub=100000, pp=False)
        print(rstr)
        if res:
            res, rstr = trial_division(p, ub=max(2000000, args.t))
            print(rstr)
            if res:
                res, rstr = proth_test(args.k, args.n, 20)
                print(rstr)
                if not res:
                    if args.f:
                        print("Warn: this might take a while")
                        res, rstr = p_minus_1_factorization(p, B=args.B)
                        print(rstr)
if __name__ == "__main__":
    main()
