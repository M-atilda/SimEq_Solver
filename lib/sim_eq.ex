defmodule SimEq do
  @moduledoc """
  this module provides kinds of methods
    to solve simultinious equations
  """

  def get_error(matrix, inhom_vector, result) do
    result_matrix = Matrix.generate(for val <- result do [val] end)
    %Matrix{contents: c} = Matrix.multiply(matrix, result_matrix)
    Enum.zip(c, inhom_vector)
    |> Enum.map(fn {[f|_], s} -> f-s end)
  end

  def get_error_norm(matrix, inhom_vector, result) do
    get_error(matrix, inhom_vector, result)
    |> Enum.reduce(0, fn (val, acm) -> (val*val)+acm end)
    |> :math.sqrt
  end

end
