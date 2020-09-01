function [val,varargout] = makeBruteForceStateTrajectories(ctrlIp,nSteps,x0,lb,ub)

% function that determines every possible state trajectory given initial
% state, discretized control, number of steps, and state bounds

% create all possible combination of control trajectories
uTraj = makeBruteForceControlCombinations(ctrlIp,nSteps);
% number of trajectories
nComb = size(uTraj,1);
% calculate determine all mean elevation trajectories
val = nan(nComb,nSteps);
val(:,1) = x0 + uTraj(:,1);
for ii = 1:nComb
    for jj = 2:nSteps
        val(ii,jj) = val(ii,jj-1) + uTraj(ii,jj);
    end
end
% eliminate trajectories going outside bounds
belowBounds = val<lb;
aboveBounds = val>ub;
outsideBounds = or(belowBounds,aboveBounds);
val(sum(outsideBounds,2)>0,:) = [];
uTraj(sum(outsideBounds,2)>0,:) = [];
varargout{1} = uTraj;

end