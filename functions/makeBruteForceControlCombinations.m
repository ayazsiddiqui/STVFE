function val = makeBruteForceControlCombinations(ctrlIp,nSteps)

% number of elements in control vector
nu = numel(ctrlIp);
% preallocate
val = NaN(nu^nSteps,nSteps);
% initialize
ii = 1;
jj = nSteps;
kk = 1;
% this weird loop that somehow works
while jj >= 1
    while ii <= nu^nSteps
        while kk <= nu
            is = nu^(nSteps-jj);
            val(ii:ii+is-1,jj) = ctrlIp(kk);
            ii = ii+is;
            kk = kk+1;
        end
        kk=1;
    end
    ii = 1;
    jj = jj-1;
end

end