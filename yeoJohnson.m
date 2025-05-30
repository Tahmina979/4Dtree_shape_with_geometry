function [transformed_data, lambda] = yeoJohnson(data, lambda)
  % Yeo-Johnson Transformation
  %
  % Args:
  %   data: The input data vector. Can contain positive and negative values.
  %   lambda: The transformation parameter. If not provided, it is estimated.
  %
  % Returns:
  %   transformed_data: The transformed data vector.
  %   lambda: The used transformation parameter.


  if nargin < 2 || isempty(lambda)
    % Estimate lambda using maximum likelihood estimation (MLE) - this part requires optimization
    %  A simple approach (not necessarily the most efficient MLE) involves a grid search.
    lambda_range = -2:0.1:2;
    log_likelihoods = zeros(size(lambda_range));
    for i = 1:length(lambda_range)
      log_likelihoods(i) = sum(log(yeoJohnson_func(data, lambda_range(i))));
    end
    [~, idx] = max(log_likelihoods);
    lambda = lambda_range(idx);

  end

  transformed_data = yeoJohnson_func(data, lambda);


end


function transformed_values = yeoJohnson_func(x, lambda)
  % Helper function to apply the Yeo-Johnson transformation for a given lambda.
  transformed_values = zeros(size(x));
  
  pos_idx = x>=0;
  neg_idx = ~pos_idx;

  if lambda == 0
    transformed_values(pos_idx) = log(x(pos_idx) + 1);
    transformed_values(neg_idx) = -log(-x(neg_idx) + 1);
  else
    transformed_values(pos_idx) = ((x(pos_idx) + 1).^lambda - 1) / lambda;
    transformed_values(neg_idx) = -( (1 - x(neg_idx)).^(2-lambda) -1) / (2-lambda);
  end

end
