% function [L,D] = cholesky(A)
function A = cholesky(A)
% in-place cholesky decomposition of A = LDL'. A must be symmetric and
% positive-definite. After calculations A contains D on the diagonal and L
% in its lower triangular part -- the diagonal entries of L are 1.

% L = diag(ones(size(A,1),1));
% D = zeros(size(A,1),1);
% for j = 1:size(A,1)
%     D(j) = A(j,j);
%     for k = 1:j-1
%         D(j) = D(j) - L(j,k)^2*D(k);
%     end
%     
%     for i = j+1:size(A,1)
%         L(i,j) = A(i,j);
%         for k = 1:j-1
%             L(i,j) = L(i,j) - L(i,k)*L(j,k)*D(k);
%         end
%         L(i,j) = L(i,j)/D(j);
%     end
% end

for j = 1:size(A,1)
    for k = 1:j-1
        A(j,j) = A(j,j) - A(j,k)^2*A(k,k);
    end
    
    for i = j+1:size(A,1)
        for k = 1:j-1
            A(i,j) = A(i,j) - A(i,k)*A(j,k)*A(k,k);
        end
        A(i,j) = A(i,j)/A(j,j);
    end
end