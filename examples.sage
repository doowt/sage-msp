from sage.all import *
load("cas.sage")
load("msp.sage")

print("Examples")

print("\nQ_2 Access Structure CAS_1: gamma_minus = {{0, 1}, {0, 2}, {0, 3}, {1, 2, 3}}\n")
cas_1 = CAS(matrix(QQ, [[1, 1, 0, 0], [1, 0, 1, 0], [1, 0, 0, 1], [0, 1, 1, 1]]), 0)
cas_1.print_all()

print("\n\n\nReplicated Secret-Sharing on CAS_1\n")
replicated_msp = Replicated_MSP(cas_1, QQ)
replicated_msp.print_all()

print("\n\n\nMultiplicative MSP from CAS_1\n")
mult_msp = Mult_MSP(replicated_msp)
mult_msp.print_all()

print("\nAccess Structure of product:\n")
mult_msp.mult_AS.print_all()


print("\n\n\nThreshold Access Structure CAS_2: n = 5, t = 2\n")
cas_2 = Threshold_CAS(5,2)
cas_2.print_all()

print("\n\n\nShamir's Secret-Sharing on CAS_2\n")
shamir_msp = Shamir_MSP(cas_2, QQ)
shamir_msp.print_all()

print("\n\n\nCompute Access Structure from MSP\n")
test_msp = MSP(CAS(matrix(QQ, [[1, 1, 1, 1]]), 0), QQ)
test_msp.M = matrix(QQ, [[1,  0,  1,  1,  1,  1],
             		 [1,  1,  0,  1,  1,  1],
             		 [1,  1,  1,  0,  1,  1],
			 [0,  1,  0,  0,  0,  0],
			 [1,  1,  1,  1,  0,  0],
			 [0,  0,  1,  0,  0,  0],
			 [0,  0,  0,  0,  1,  0],
			 [0,  0,  0,  1,  0,  0],
                         [0,  0,  0,  0,  0,  1]])
test_msp.rowmap = matrix(QQ, [[0], [0], [0], [1], [1], [2], [2], [3], [3]])
test_msp.target = matrix(QQ, [1] * 6)
test_msp.compute_cas()
test_msp.compute_all_msp_properties()
test_msp.AS.print_all()
print()
test_msp.print_all()

print("\n\n\nTransform MSP to Multiplicative MSP\n")
test_mmsp = Mult_MSP(test_msp)
test_mmsp.print_all()

print("\nAccess Structure of product:\n")
test_mmsp.mult_AS.print_all()
