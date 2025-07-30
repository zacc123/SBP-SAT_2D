"""
Testing cases for 2D SBP Ops
"""

include("./diagonal_sbp.jl")
include("./sbp_ops.jl")

include("./helper_functions.jl")

"""
Tests to confirm behavior of SBP Operators matches 
Erickson + Dunham (2014) for constant μ
Almquist + Dunahm (2021) for var μ
"""
function op_checks()

    # const ops chec
    
    # define basic grid
    R0, RN, DR = (-1, 1, 0.1)
    S0, SN, DS = (-1, 1, 0.1)

    NR = Int((RN - R0) / DR - 1)
    NS = Int((SN - S0) / DS - 1)

    # set order of operator 
    p = 2
    (D2s, D2r, Is, Ir, Hs, Hr, HIs, HIr, BSs, BSr, mu, Ef, Er, Es, Ed) = sbp_operators(p, S0, SN, R0, RN, NS, NR, DS, DR)


    # test vector:
    stacked_x = zeros((NR + 1) * (NS + 1))
    set_vector!(stacked_x, NR, NS)
    
    print("\t\tTesting Lift Operators....\n")
    test_vector_lift(NR, NS, stacked_x, Ef, Er, Es, Ed)
    print("\t\t[PASS] Test Lift Operator\n")

end


"""
Test that Ef, etc. Operators lift out the correct data
    Assumes l and r side are stacking directions (S)
            t and b are (R)
"""
function test_vector_lift(NR, NS, stacked_x, lift_t, lift_b, lift_l, lift_r)

    test_t = zeros((NR + 1) * (NS + 1))
    test_b = zeros((NR + 1) * (NS + 1))
    test_l = zeros((NR + 1) * (NS + 1))
    test_r = zeros((NR + 1) * (NS + 1))

    lift_r_vectors!(stacked_x, test_t, test_b, NR, NS, 1, NS+1)
    lift_s_vectors!(stacked_x, test_l, test_r, NR, NS, 1, NR+1)

    # print(test_t, "\n\n", lift_t * stacked_x, "\n\n")
    # print(test_b, "\n\n", lift_b * stacked_x, "\n\n")
    # print(test_l, "\n\n", lift_l * stacked_x, "\n\n")
    # print(test_r, "\n\n", lift_r * stacked_x, "\n\n")

    @assert test_t == lift_t * stacked_x
    @assert test_b == lift_b * stacked_x
    @assert test_l == lift_l * stacked_x
    @assert test_r == lift_r * stacked_x

    return nothing
end

"""
Initializes a vector so that R goes 1 -> NR + 1, stacked on S
  ^   
S |
R  -  >

    [
     [1, 2, 3, 4, 5],
     [1, 2, 3, 4, 5],
     [1, 2, 3, 4, 5],
     [1, 2, 3, 4, 5]
     ]
"""
function set_vector!(c, NR, NS)
    
    for i in 1:NS+1

        for j in 1:NR+1

            c[(i - 1)*(NR+1) + j] = j
        
        end

    end

end


############# MAIN Function to Run Tests ################
# bear with the tabbing plzzzzz
function main()
    print("\nRunning Test Suite....\n")

        print("\tRunning Basic Operator Checks:\n")
            op_checks()
        print("\t[PASS] Basic Operator Checks:\n")

    print("[PASS] Test Suite\n")
end

main()