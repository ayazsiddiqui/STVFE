function val = makeBruteForceCombinations(ctrlIp,nSteps)

% function that determines every possible control sequence that can be
% executed given a control input vector and number of steps

nu = numel(ctrlIp);

val = NaN(nu^nSteps,nSteps);

ii = 1;
jj = nSteps;
kk = 1;

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