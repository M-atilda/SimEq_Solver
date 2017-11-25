defmodule IterationTest do
  use ExUnit.Case
  doctest SimEq

  @sample_m Matrix.generate([[2.0, -1.0, 0.0],
                             [-1.0, 2.0, -1.0],
                             [0.0, -1.0, 2.0]])
  @sample_b [1.0, 0.0, 1.0]
  @sample_x [0.0, 0.0, 0.0]

  test "calc with jacob method" do
    {ite, result} = SimEq.Iteration.solve_jacob(
      @sample_m,
      @sample_b,
      @sample_x)

    IO.puts ""
    IO.puts "iteration: #{ite}"
    IO.puts "result:"
    Enum.map(result, fn val -> IO.puts "  #{val}" end)
  end

  test "calc with gauss seidel method" do
    {ite, result} = SimEq.Iteration.solve_gauss_seidel(
      @sample_m,
      @sample_b,
      @sample_x)

    IO.puts ""
    IO.puts "iteration: #{ite}"
    IO.puts "result:"
    Enum.map(result, fn val -> IO.puts "  #{val}" end)
  end

end
