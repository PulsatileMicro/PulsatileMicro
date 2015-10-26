function stop = pswIteraRes(optimValues,state)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

stop = false; % This function does not stop the solver
switch state
	case 'init'
        timeNow=datestr(now);
        fid = fopen('AdapCoeff_Wall_res.txt','wt');
        fprintf(fid, timeNow);
        fid = fclose(fid);
        record = optimValues.swarm;
        save('AdapCoeff_Wall_res.txt','record','-ascii','-v6','-append');
        tic;
    case 'iter'
        elaseTime=toc;
        disp(['particleswarms第',num2str(optimValues.iteration),'次迭代运行时间',num2str(elaseTime)]);
        record=[elaseTime,optimValues.bestfval,optimValues.funccount];
        iteraCoeff = optimValues.bestx;
        save('AdapCoeff_Wall_res.txt','record','iteraCoeff','-ascii','-v6','-append');
        tic;
    case 'done'
        % No cleanup necessary
    case 'interrupt'
        % TODO
        stop = true;
end