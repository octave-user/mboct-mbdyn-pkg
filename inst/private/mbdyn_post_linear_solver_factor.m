function x = mbdyn_post_linear_solver_factor(A, b, options)
  x = fem_sol_factor(A, options) \ b;
endfunction
