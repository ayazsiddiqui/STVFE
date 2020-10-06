classdef plottingSuite
    %PLOTTINGSUITE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    
    
    %% path related plots
    methods
        % plot dome given radius
        function val = plotDome(~,radius)
            % get constants
            r = radius;
            linWidth = 0.5;
            lnType = ':';
            grayRGB = 128/255.*[1 1 1];
            % make longitude and latitude fine grids
            longFine = -90:1:90;
            latFine = -0:1:90;
            % make longitude and latitude coarse grids
            stepSize = 30;
            longCoarse = longFine(1):stepSize:longFine(end);
            latCoarse = latFine(1):stepSize:latFine(end);
            val = gobjects;
            kk = 1;
            % plot longitude lines
            for ii = 1:numel(longCoarse)
                X = r*cosd(longCoarse(ii)).*cosd(latFine);
                Y = r*sind(longCoarse(ii)).*cosd(latFine);
                Z = r*sind(latFine);
                val(kk) = plot3(X,Y,Z,lnType,'linewidth',linWidth,'color',grayRGB);
                hold on;
                kk = kk+1;
            end
            % plot latitude lines
            for ii = 1:numel(latCoarse)
                X = r*cosd(longFine).*cosd(latCoarse(ii));
                Y = r*sind(longFine).*cosd(latCoarse(ii));
                Z = r*sind(latCoarse(ii))*ones(size(longFine));
                val(kk) = plot3(X,Y,Z,lnType,'linewidth',linWidth,'color',grayRGB);
                kk = kk+1;
            end
            % view angle
            view(100,35);
        end
        
        % plot path given path parameters
        function val = plotPath(~,radius,pathBasis,prevPath)
            % local variables
            s =linspace(0,2*pi,201);
            % path basis parameters
            pC = getPathCoords(pathBasis(1),...
                pathBasis(2),pathBasis(3),radius,s);
            % get lemniscate coordinates
            lemX = pC(1,:);
            lemY = pC(2,:);
            lemZ = pC(3,:);
            % plot path
            if nargin == 3
                val = plot3(lemX,lemY,lemZ,'k-',...
                    'linewidth',1);
            else
                prevPath.XData = lemX;
                prevPath.YData = lemY;
                prevPath.ZData = lemZ;
                
            end
            view(100,35);
        end
        
    end
    
    
    %% vector related plotting methods
    % make a quiver plot for vector
    methods
        function val = plotQuiver(~,location,vecVal,qv)
            vecVal = vecVal./norm(vecVal);
            if nargin == 3
                val = quiver3(location(1),location(2),location(3),...
                    vecVal(1),vecVal(2),vecVal(3),...
                    'MaxHeadSize',0.6,...
                    'linewidth',0.8);
            else
                val = qv;
                val.XData = location(1);
                val.YData = location(2);
                val.ZData = location(3);
                val.UData = vecVal(1);
                val.VData = vecVal(2);
                val.WData = vecVal(3);
            end
        end
    end
    
    
    %% plot simulation vectors
    methods
        % plot an axis system
        function val = plotAxes(obj,location,vecVal,prevAxis)
            val = gobjects([3,1]);
            rgbColors = 1/255*[228,26,28        % red
                77,175,74                       % green
                55,126,184];                    % blue
            if nargin == 3
                for ii = 1:3
                    val(ii) = obj.plotQuiver(location,vecVal(:,ii));
                    val(ii).Color = rgbColors(ii,:);
                    hold on;
                end
            else
                for ii = 1:3
                    val = prevAxis;
                    val(ii) = ...
                        obj.plotQuiver(location,vecVal(:,ii),prevAxis(ii));
                end
            end
            
        end
        
    end
    
    
    %% animation methods
    methods
        % make animation
        function val = makeAnimation(obj,G_rCM,G_bdy,pathBasis,varargin)
            time = G_rCM.Time;
            nPoints = numel(time);
            % parse input
            pp = inputParser;
            addParameter(pp,'waitforbutton',true,@islogical);
            addParameter(pp,'G_vCM',timeseries(zeros(3,nPoints),time));
            parse(pp,varargin{:});
            dRadius = vecnorm(squeeze(G_rCM.Data(:,1,:)));
            % plot the dome
            pDome = obj.plotDome(dRadius(1));
            % annotations
            grid on; hold on;
            xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
            % plopath
            pPath = obj.plotPath(dRadius(1),pathBasis.Data(:,1,1));
            % plot position vec
            pPos = plot3([0 G_rCM.Data(1,1,1)],[0 G_rCM.Data(2,1,1)],...
                [0 G_rCM.Data(3,1,1)],'b-');
            % plot body axes
            pBdy = obj.plotAxes(G_rCM.Data(:,:,1),G_bdy.Data(:,:,1));
            for ii = 1:3
                pBdy(ii).AutoScaleFactor = dRadius(1)*0.15;
            end
            % plot velocity
            pVel = obj.plotQuiver(G_rCM.Data(:,:,1),...
                pp.Results.G_vCM.Data(:,:,1));
            pVel.Color = [255,127,0]/255;
            pVel.AutoScaleFactor = dRadius(1)*0.15;
            % title
            titleTxt = sprintf('Time = %.2f sec', G_rCM.Time(ii));
            title(titleTxt);
            % video
            val = VideoWriter('animVid.avi','Uncompressed AVI');
            val.FrameRate = 1/(G_rCM.Time(2) - G_rCM.Time(1));
            open(val);
            frame = getframe(gcf);
            writeVideo(val,frame);
            
            for ii = 2:nPoints
                if numel(unique(round(dRadius))) ~= 1
                    delete(pDome);
                    pDome = obj.plotDome(dRadius(ii));
                end
                if numel(unique(squeeze(pathBasis.Data)','rows')) ~= 3
                    delete(pPath);
                    pPath = obj.plotPath(dRadius(ii),...
                        pathBasis.Data(:,1,ii),pPath);
                end
                % update position
                pPos.XData = [0 G_rCM.Data(1,1,ii)];
                pPos.YData = [0 G_rCM.Data(2,1,ii)];
                pPos.ZData = [0 G_rCM.Data(3,1,ii)];
                % body axes
                pBdy = obj.plotAxes(G_rCM.Data(:,1,ii),...
                    G_bdy.Data(:,:,ii),pBdy);
                % velocity
                pVel = obj.plotQuiver(G_rCM.Data(:,:,ii),...
                pp.Results.G_vCM.Data(:,:,ii),pVel);
                % title
                titleTxt = sprintf('Time = %.2f sec', time(ii));
                title(titleTxt);
                
                if pp.Results.waitforbutton
                    waitforbuttonpress;
                end
                
                frame = getframe(gcf);
                writeVideo(val,frame);
                
            end
            
            close(val); % Close the movie file
            
        end
    end
end


