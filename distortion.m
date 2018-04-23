function error=distortion(original,reconstructed)

%estimation of distortion error between the original data and the
%reconstructed via the nueral gas algorithm
%plsease see equation 2

[sensor time]=size(original); % sensors x time

terror=zeros(1,sensor);

for s=1:sensor % sensors
    numerator=0;
    denominator=0;

    for t=1:time %time
        numerator = numerator + (original(sensor,time) - reconstructed(sensor,time))^2;
        denominator = denominator + (original(sensor,time) - mean(original(sensor,:)))^2;
    end
    terror(s)=numerator/denominator;
end

error=mean(terror);

