% N = length_MSD; Note to self
% epsilon = loc_uncert; Note to self
% function sigma = myMichaletFun(N, epsilon_val, alpha_val)
function MSD_data = myMichaletFun(MSD_data)
for i = 1:length(MSD_data)
%From the paper, the varience of the mean squared displacement is:
    N = length(MSD_data(i).MSD);
    alpha_val = MSD_data(i).alpha_val;
    epsilon_val = MSD_data(i).epsilon_val;
    sigma = zeros(1,N);
    for n = 1:N
        K = N - n; 
        if n <= K
            sigma(n) = sqrt(((n/(6*(K^2)))*(((4*(n^2)*K)+(2*K)-(n^3)+n)*(alpha_val^2))) + ((1/K)*((2*n*alpha_val*epsilon_val)+(1+(0.5*(1-(n/K)))*epsilon_val^2))));
        elseif n > K
            sigma(n) = sqrt(((1/(6*K))*(((6*(n^2)*K)-(4*n*(K^2))+(4*n)+(K^3)- K)*(alpha_val^2)))+ ((1/K)*((2*n*alpha_val*epsilon_val)+(epsilon_val^2))));
        end
    end
    MSD_data(i).MSD_MichaeletError = sigma;
    disp(['Michaeletfun#', num2str(i)]);
end
end

%     sigma = zeros(1,N);
%     for n = 1:N
%         K = N - n; 
%         if n <= K
%             sigma(n) = sqrt(((n/(6*(K^2)))*(((4*(n^2)*K)+(2*K)-(n^3)+n)*(MSD_data(i).alpha_val(i)^2))) + ((1/K)*((2*n*MSD_data(i).alpha_val(i)*MSD_data(i).epsilon_val(i))+(1+(0.5*(1-(n/K)))*MSD_data(i).epsilon_val(i)^2))));
%         elseif n > K
%             sigma(n) = sqrt(((1/(6*K))*(((6*(n^2)*K)-(4*n*(K^2))+(4*n)+(K^3)- K)*(MSD_data(i).alpha_val(i)^2)))+ ((1/K)*((2*n*MSD_data(i).alpha_val(i)*MSD_data(i).epsilon_val(i))+(MSD_data(i).epsilon_val(i)^2))));
%         end
%     end
%     MSD_data(i).MSD_MichaeletError = sigma;
%     disp(['Michaeletfun#', num2str(i)]);


%     N = length(MSD_data(i).MSD);
% From David's Report
% sigma = zeros(1,N);
% for n = 1:N
%     K = N - n; 
%     if n <= K
%         sigma(n) = sqrt(((n/(6*(K^2)))*(((4*(n^2)*K)+(2*K)-(n^3)+n)*(alpha_val^2))) + ((1/K)*((2*n*alpha_val*epsilon_val)+(1+(0.5*(1-(n/K)))*epsilon_val^2))));
%     elseif n > K
%         sigma(n) = sqrt(((1/(6*K))*(((6*(n^2)*K)-(4*n*(K^2))+(4*n)+(K^3)- K)*(alpha_val^2)))+ ((1/K)*((2*n*alpha_val*epsilon_val)+(epsilon_val^2))));
%     end
% end
% end
% end