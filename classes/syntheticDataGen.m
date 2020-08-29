classdef syntheticDataGen
    %% properties
    properties
        spatialKernel
        temporalKernel
    end
    
    properties (SetAccess = private)
        temporalCovariance
        spatialCovariance
    end
    
    %% methods
    methods
        % % % %         constructor
        function obj = syntheticDataGen(varargin)
            % % parse input
            pp = inputParser;
            addParameter(pp,'spatialKernel','squaredExponential',@ischar);
            addParameter(pp,'temporalKernel','exponential',@ischar);
            parse(pp,varargin{:});
            spKer = pp.Results.spatialKernel;
            tKer = pp.Results.temporalKernel;
            % % list of programmed kernels
            progKernels = {'exponential','squaredExponential'};
            % % check if the selected kernels are programmed
            % reference: http://www.gaussianprocess.org/gpml/chapters/RW4.pdf
            % Eqns. (4.9) and (4.18)
            % % temporal kernel
            obj.temporalKernel = tKer;
            switch obj.temporalKernel
                case 'exponential'
                    obj.temporalCovariance = @(t1,t2,timeScale)...
                        exp(-abs(t2-t1)/timeScale);
                case 'squaredExponential'
                    obj.temporalCovariance = @(t1,t2,timeScale)...
                        exp(-0.5*((t2-t1)^2)/timeScale^2);
                otherwise % error message
                    error(['Only ',repmat('"%s", ',1,numel(progKernels)),...
                        'are valid entries for temporal kernel. '...
                        'You enterd "%s". \n'],progKernels{:},tKer);
            end
            % % spatial kernel
            obj.spatialKernel = spKer;
            switch obj.spatialKernel
                case 'exponential'
                    obj.spatialCovariance = @(s1,s2,xLs,yLs,zLs)...
                        exp(-sum(abs(s1(:)-s2(:))./[xLs;yLs;zLs]));
                case 'squaredExponential'
                    obj.spatialCovariance = @(s1,s2,xLs,yLs,zLs)...
                        exp(-0.5*(s1-s2)'*(eye(3)...
                        ./[xLs;yLs;zLs].^2)*(s1-s2));
                otherwise % error message
                    error(['Only ',repmat('"%s", ',1,numel(progKernels)),...
                        'are valid entries for spatial kernel. '...
                        'You enterd "%s". \n'],progKernels{:},spKer);
            end
        end
        
        % % % %         get the colored field
        function val = getColoredData(obj,x,y,z,t,xLs,yLs,zLs,tS)
            % % number of elements in the different grids
            xNp = numel(x);
            yNp = numel(y);
            zNp = numel(z);
            tNp = numel(t);
            % % make a grid of spatial locations,
            % % didn't like meshgrid so I'm using for loops
            % total number of spatial locations
            spTotNp = xNp*yNp*zNp;
            % preallocate matrix
            spatialLocs = NaN(3,spTotNp);
            cc = 1;
            % % fill it up
            for ii = 1:xNp
                for jj = 1:yNp
                    for kk = 1:zNp
                        spatialLocs(:,cc) = [x(ii);y(jj);z(kk)];
                        cc = cc+1;
                    end
                end
            end
            % % make the upper triangular spatial covariance matrix
            spCov = NaN(spTotNp);
            for ii = 1:spTotNp
                for jj = ii:spTotNp
                    spCov(ii,jj) = obj.spatialCovariance(spatialLocs(:,ii),...
                        spatialLocs(:,jj),xLs,yLs,zLs);
                end
            end
            % % form the total covariance matrix and add small number
            spCov = spCov + triu(spCov,1)' + 0.0001*eye(spTotNp);
            % % make the upper triangular temporal covariance matrix
            tCov = NaN(tNp);
            for ii = 1:tNp
                for jj = ii:tNp
                    tCov(ii,jj) = obj.temporalCovariance(t(:,ii),...
                        t(:,jj),tS);
                end
            end
            % % form the total covariance matrix and add small number
            tCov = tCov + triu(tCov,1)' + 0.0001*eye(tNp);
            % % take the cholesky decompositions of the two matrices
            Lsp = chol(spCov);
            Lt = chol(tCov);
            % % take random number samples
            warning('Using random numbers drawn from a normal distribtion.')
            samp = randn(spTotNp,tNp);
            % % colored signal
            coloredSignal = (Lsp*(Lt*samp')');
            % % outputs
            val.coloredSignal = coloredSignal;
            val.dsgnSpace = cat(1,repmat(spatialLocs,1,1,tNp),...
                reshape((t(:).*ones(1,spTotNp))',1,spTotNp,tNp));
            
        end
        
        % % % %         visualization function
        function val = visualize(obj,coloredSignal,dsgnSpace)
            % % marker size
            markerSize = 15;
            % % figure props
            figure(1);
            set(gcf,'position',[1268 466 1.15*[560 420]]);
            F = struct('cdata',uint8(zeros(840,1680,3)),'colormap',[]);            
            
            % normalize between 0 and 1
            coloredSignalNorm = (-min(coloredSignal,[],'all')+coloredSignal)...
                ./range(coloredSignal,'all');
            for ii = 1:size(dsgnSpace,3)
                                
                if ii == 1
                    xlabel('x')
                    ylabel('y')
                    zlabel('z')
                    grid on
                    hold on
                    view(-35,30)
                else
                    delete(pq);
                end
                pq = scatter3(dsgnSpace(1,:,ii),dsgnSpace(2,:,ii),...
                    dsgnSpace(3,:,ii),markerSize,coloredSignalNorm(:,ii)',...
                    'filled');
                co = colorbar;
                co.Label.String = 'Normalized flow speed';
                co.Label.Position(1) = -1;

                % % title
                txt = sprintf('Normalized time = %.2f',...
                    dsgnSpace(4,1,ii)/dsgnSpace(4,1,end));
                title(txt);
                
                ff = getframe(gcf);
                F(ii).cdata = ff.cdata;
            end
            % % output the frame
            val = F;
            
        end
        
    end
end
