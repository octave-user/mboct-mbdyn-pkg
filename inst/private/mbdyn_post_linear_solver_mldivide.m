function x = mbdyn_post_linear_solver_mldivide(A, b, options)
  x = full(A \ b);
endfunction
