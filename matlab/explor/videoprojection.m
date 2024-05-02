classdef videoprojection < projection
    % function P = videoprojection(G,videoObj[,nframe])
    %---
    % Project a video object representing a 4D array (x, y, color, time)
    % into 3D images (x, y, color)
    
    % communication properties
    properties (SetAccess='private')
        videoObj
        rgb
        button
    end
            
    % Constructor, destructor, and events
    methods
        function P = videoprojection(G,videoObj,nframe)
            fn4D_dbstack
            P = P@projection(G,[1 2]);
            % check video object
            if nargin<3
                disp('read frames to obtain video''s number of frames')
                nframe = 0;
                step = 10000;
                while true
                    k = nframe + step;
                    try
                        videoObj.read(k);
                        nframe = k;
                    catch
                        if step == 1
                            break
                        else
                            step = step/10;
                        end
                    end
                end
                disp(['number of frames: ' num2str(nframe)])
            end       
            frame = videoObj.read(1);
            P.videoObj = videoObj;
            P.rgb = (size(frame,3)==3);
            % set dimsplus
            if P.rgb
                P.dimsplus = 3;
            end
            % set data
            if P.rgb
                P.G.sizes(1:4) = [videoObj.Width videoObj.Height 3 nframe];
            else
                P.G.sizes(1:3) = [videoObj.Width videoObj.Height nframe];
            end
            P.videoObj = videoObj;   
            P.sliceframe()
        end
    end
        
    % Update upon events
    methods
        function updateDown(P,~,evnt)
            fn4D_dbstack(['G2S ' evnt.flag])
            g=P.G;
            % not following selections in the orthogonal dimensions, but
            % only 'ijkl'!
            switch evnt.flag
                case 'ijkl'
                    % update frame only if change in dimension 3+P.rgb
                    ijklold = evnt.oldvalue;
                    ijklchg = (g.ijkl~=ijklold);
                    if ijklchg(3+P.rgb)
                        P.sliceframe()
                    end
                    return
                case 'selection'
                    if all(ismember(evnt.dims,P.datadims.nodisplay)) ...
                             || isempty(evnt.dims)
                         return
                    end
            end
            updateDown@projection(P,[],evnt)
        end
    end
    
    % Frame slicing
    methods (Access='private')
        function sliceframe(P)
            global kk
            fn4D_dbstack
            kframe = P.G.ijkl(3+P.rgb);
            kk = [kk kframe];
            frame = permute(P.videoObj.read(kframe),[2 1 3]);
            P.slice = struct('active',true,'data',frame);
        end
    end
    
    % Local menu
    methods
        function initlocalmenu(P,button)
            P.button = button;
        end
    end
end
