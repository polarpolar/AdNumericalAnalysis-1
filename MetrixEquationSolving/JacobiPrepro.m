function [ M, b ] = JacobiPrepro( M, b )
%JACOBIPREPRO Summary of this function goes here
%   This is the Jacobi pre-process approach to optimize the input metrix by
%   reducing the condition number.

    M = diag(diag(M)) \ M;
    b = diag(diag(M)) \ b;


end

