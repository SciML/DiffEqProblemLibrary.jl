module BVProblemLibrary

using DiffEqBase, Markdown

include("linear.jl")
include("nonlinear.jl")

# Linear BVP Example Problems
export  prob_bvp_linear_1, prob_bvp_linear_2, prob_bvp_linear_3, prob_bvp_linear_4, prob_bvp_linear_5,
        prob_bvp_linear_6, prob_bvp_linear_7, prob_bvp_linear_8, prob_bvp_linear_9, prob_bvp_linear_10,
        prob_bvp_linear_11, prob_bvp_linear_12, prob_bvp_linear_13, prob_bvp_linear_14, prob_bvp_linear_15,
        prob_bvp_linear_16, prob_bvp_linear_17, prob_bvp_linear_18

# Nonlinear BVP Example Problems
export  prob_bvp_nonlinear_1, prob_bvp_nonlinear_2, prob_bvp_nonlinear_3, prob_bvp_nonlinear_4, prob_bvp_nonlinear_5,
        prob_bvp_nonlinear_6, prob_bvp_nonlinear_7, prob_bvp_nonlinear_8, prob_bvp_nonlinear_9, prob_bvp_nonlinear_10,
        prob_bvp_nonlinear_11, prob_bvp_nonlinear_12, prob_bvp_nonlinear_13, prob_bvp_nonlinear_14, prob_bvp_nonlinear_15


end # module BVProblemLibrary
