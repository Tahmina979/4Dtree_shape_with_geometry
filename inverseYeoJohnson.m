function original_data = inverseYeoJohnson(transformed_data, lambda)
  % Inverse Yeo-Johnson Transformation
  %
  % Args:
  %   transformed_data: The transformed data vector.
  %   lambda: The transformation parameter used in the forward transformation.
  %
  % Returns:
  %   original_data: The original data vector.

  original_data = zeros(size(transformed_data));
  pos_idx = transformed_data >= 0;
  neg_idx = ~pos_idx;

  if lambda == 0
    original_data(pos_idx) = exp(transformed_data(pos_idx)) - 1;
    original_data(neg_idx) = -exp(-transformed_data(neg_idx)) + 1;
  else
    original_data(pos_idx) = (lambda * transformed_data(pos_idx) + 1).^(1/lambda) -1;
    original_data(neg_idx) = 1 - ( (2 - lambda) * (-transformed_data(neg_idx)) + 1).^((1/(2-lambda)));
  end

end