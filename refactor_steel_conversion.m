function steelParam = refactor_steel_conversion(param)
    % Steel model conversion.
    %   Calculate totalAs, steelratio
    
    % load the parameters
    handles = param.handles;
    name = param.name;
    s = param.s;
    ss = param.ss;
    b = param.b;
    h = param.h;
    cover = param.cover;
    
%     p=1;
%     for i=1:1:mat(1)
%         areaAs(i)=0;j=1;
%         for j=1:1:(mat(2)-2)
%             switch ss(i,j)
%                 case 0
%                     As=0;
%                 case 3
%                     As=0.71;
%                 case 4
%                     As=1.27;
%                 case 5
%                     As=1.99;
%                 case 6
%                     As=2.87;
%                 case 7
%                     As=3.87;
%                 case 8
%                     As=5.07;
%                 case 9
%                     As=6.47;
%                 case 10
%                     As=8.14;
%                 case 11
%                     As=10.07;
%                 case 12
%                     As=12.19;
%                 case 14
%                     As=14.52;
%                 case 16
%                     As=19.79;
%                 case 18
%                     As=25.79;
%             end
%             sAs(i,j)=As;
%             if As>0
%                 ssAs(1,p)=As;
%                 p=p+1;
%             end
%             areaAs(i)=areaAs(i)+sAs(i,j);
%             As=0;
%         end
%     end
%     steelnum=(p-1);
%     totalAs=0;dtmax=1;
%     for n=1:1:mat(1)
%         d(n)=s(n,1);
%         totalAs=totalAs+areaAs(n);
%         if d(n)>dtmax
%             dtmax=d(n);
%         end
%     end
    As = getSteelModel();
    ssAs = reshape(As(ss + 1)', 1, []);    % plus one is because matlab indexed from one instead of zero
    ssAs(ssAs == 0) = [];
    steelnum = length(ssAs);
    totalAs = sum(ssAs);
    dtmax = max(s(:,1));

    steelratio = roundn(totalAs / (b * h), -4);
    set(handles.text28, 'String', name);
    set(handles.text4, 'String', cover);
    set(handles.text6, 'String', steelratio);
    
    steelParam.handles = handles;
    steelParam.ssAs = ssAs;
    steelParam.steelnum = steelnum;
    steelParam.totalAs = totalAs;
    steelParam.dtmax = dtmax;
end

function As = getSteelModel()
    % indexed from 1 ~ 19 presenting model 0 ~ 18
    As = [0, ...
        nan, nan, 0.71, 1.27, 1.99, 2.87, 3.87, 5.07, 6.47, 8.14, ...
        10.07, 12.19, nan, 14.52, nan, 19.79, nan, 25.79];
end