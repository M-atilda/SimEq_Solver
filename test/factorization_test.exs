defmodule FactorizationTest do
  use ExUnit.Case
  doctest SimEq

  @sample_m SimEq.Matrix.generate([[2.0, 3.0, -1.0],
                                   [4.0, 4.0, -3.0],
                                   [2.0, -3.0, 1.0]])
  @sample_b [5.0, 3.0, -1.0]

  # test "calc with gaussian method" do
  #   assert SimEq.Gauss.solve_gaussian_eliminate(@sample_m, @sample_b, true, true) == [1.0, 2.0, 3.0]
  # end

  test "calc with LU factorization" do
    assert SimEq.Factorization.solve_factorize(@sample_m, [@sample_b])
  end

end
