function x = model_response(celltype,input,srf,ada,norm)
%
%

if norm
    x = ((sum(input(:).*srf(:)))/(sum(input(:).*ada(:))));
else
    x = (sum(input(:).*srf(:)));
end

switch celltype
    
    case 'ON'
        
        x(x < 0) = 0;
        x = abs(x);
        
    case 'OFF'
        
        x(x > 0) = 0;
        x = abs(x);
        
end