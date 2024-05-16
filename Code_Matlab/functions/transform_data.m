function [x22,mut,sdt]=transform_data(x2,DEMEAN)
% =========================================================================
% DESCRIPTION
% This function transforms a given set of series based upon the input
% variable DEMEAN. The following transformations are possible:
%
%   1) No transformation.
%   
%   2) Each series is demeaned only (i.e. each series is rescaled to have a
%   mean of 0).
%   
%   3) Each series is demeaned and standardized (i.e. each series is
%   rescaled to have a mean of 0 and a standard deviation of 1).
%   
%   4) Each series is recursively demeaned and then standardized. For a
%   given series x(t), where t=1,...,T, the recursively demeaned series
%   x'(t) is calculated as x'(t) = x(t) - mean(x(1:t)). After the
%   recursively demeaned series x'(t) is calculated, it is standardized by
%   dividing x'(t) by the standard deviation of the original series x. Note
%   that this transformation does not rescale the original series to have a
%   specified mean or standard deviation.
%
% -------------------------------------------------------------------------
% INPUTS
%           x2      = set of series to be transformed (one series per
%                     column); no missing values;
%           DEMEAN  = an integer indicating the type of transformation
%                     performed on each series in x2; it can take on the
%                     following values:
%                           0 (no transformation)
%                           1 (demean only)
%                           2 (demean and standardize)
%                           3 (recursively demean and then standardize) 
%
% OUTPUTS
%           x22     = transformed dataset
%           mut     = matrix containing the values subtracted from x2
%                     during the transformation
%           sdt     = matrix containing the values that x2 was divided by
%                     during the transformation
%
% =========================================================================
% FUNCTION

% Number of observations in each series (i.e. number of rows in x2)
T=size(x2,1);

% Number of series (i.e. number of columns in x2)
N=size(x2,2);

% Perform transformation based on type determined by 'DEMEAN'
switch DEMEAN
    
    % ---------------------------------------------------------------------
    % No transformation
    case 0
        mut=repmat(zeros(1,N),T,1);
        sdt=repmat(ones(1,N),T,1);
        x22=x2;
        
    % ---------------------------------------------------------------------
    % Each series is demeaned only
    case 1
        mut=repmat(mean(x2),T,1);
        sdt=repmat(ones(1,N),T,1);
        x22=x2-mut;
        
    % ---------------------------------------------------------------------
    % Each series is demeaned and standardized 
    case 2
        mut=repmat(mean(x2),T,1);
        sdt=repmat(std(x2),T,1);
        x22=(x2-mut)./sdt;
        
    % ---------------------------------------------------------------------
    % Each series is recursively demeaned and then standardized
    case 3
        mut=NaN(size(x2));
        for t=1:T
            mut(t,:)=mean(x2(1:t,:),1);
        end
        sdt=repmat(std(x2),T,1);
        x22=(x2-mut)./sdt; 
end





