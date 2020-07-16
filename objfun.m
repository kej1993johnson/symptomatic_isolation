function [err] = objfun(theta, t, yobs,disease, distrib)
switch disease
    case 'flu'
        switch distrib
            case 'gamma'
                err = ((gampdf(t, theta(1), theta(2)) - yobs));

             case 'normal'
                err = ((normpdf(t, theta(1), theta(2)) - yobs));
        end
        
      case 'SARS'
        switch distrib
            case 'gamma'
                err = ((gampdf(t, theta(1), theta(2)) - yobs));

             case 'normal'
                err = ((normpdf(t, theta(1), theta(2)) - yobs));
            end
              
    end



end