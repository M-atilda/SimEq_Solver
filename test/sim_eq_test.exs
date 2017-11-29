defmodule SimEqTest do
  use ExUnit.Case
  doctest SimEq

  # @tag timeout: 10800000

  # @matrix Matrix.generate([[2.0, 3.0, -1.0],
  #                          [4.0, 4.0, -3.0],
  #                          [2.0, -3.0, 1.0]])
  # @inhom [5.0, 3.0, -1.0]
  #@start [0.0, 0.0, 0.0]

  # @matrix Matrix.generate([[5, -3, -1, 0, 2, 1],
  #                          [0, -5, 2, 1, 0, 2],
  #                          [-1, 3, -6, -2, 3, 3],
  #                          [-1, 0, 3, 5, -2, -1],
  #                          [0, 3, 3, -1, 5, -4],
  #                          [2, 3, 2, 3, 2, -4]])
  # @inhom [16, 12, 12, 16, 18, 20]
  # @start [10, 10, 10, 10, 10, 10]
  # # @start [0,0,0,0,0,0]

  @matrix Matrix.generate(
    [[   2,   0,   0,   0,  -1,   0,-0.5, 0.5,   0],
     [   0,   2,   0,  -1,   0,   0, 0.5,-0.5,   0],
     [   0,   0, 2.5, 0.5,-0.5,-0.5,  -1,   0,   0],
     [   0,  -1, 0.5, 1.5,-0.5,-0.5,   0,   0,   0],
     [  -1,   0,-0.5,-0.5,   2,   0,   0,   0,-0.5],
     [   0,   0,-0.5,-0.5,   0,   2,   0,  -1, 0.5],
     [-0.5, 0.5,  -1,   0,   0,   0, 2.5,-0.5,  -1],
     [ 0.5,-0.5,   0,   0,   0,  -1,-0.5, 1.5,   0],
     [   0,   0,   0,   0,-0.5, 0.5,  -1,   0, 1.5]])
  p_3y_ex = -500 / (2.8 * :math.pow(10, 7))
  p_5y_ex = -300 / (2.8 * :math.pow(10, 7))
  #p_3y_ex = -500
  #p_5y_ex = -300
  @inhom [0,0,0,p_3y_ex,0,0,0,p_5y_ex,0]
  # @start [1,1,1,1,1,1,1,1,1]


  test "gaussian" do

    IO.puts "@gaussian"
    IO.puts "* result *"
    result = SimEq.Gauss.solve_gaussian_eliminate(@matrix, @inhom)
    Enum.map(result, fn val -> IO.puts val end)

    IO.puts "* error *"
    error = SimEq.get_error(@matrix, @inhom, result)
    Enum.map(error, fn val -> IO.puts val end)
    IO.puts "* error norm *"
    IO.puts SimEq.get_error_norm(@matrix, @inhom, result)
  end

  # test "factorization" do
  #   IO.puts "@LU factorization"
  #   IO.puts "* result *"
  #   [result|_] = SimEq.Factorization.solve_factorize(@matrix, [@inhom])
  #   Enum.map(result, fn val -> IO.puts val end)

  #   IO.puts "* error *"
  #   error = SimEq.get_error(@matrix, @inhom, result)
  #   Enum.map(error, fn val -> IO.puts val end)
  #   IO.puts "* error norm *"
  #   IO.puts SimEq.get_error_norm(@matrix, @inhom, result)
  # end


  # test "calc with jacob method" do
  #   {ite, result} = SimEq.Iteration.solve_jacob(
  #     @matrix,
  #     @inhom,
  #     @start,
  #     1000)

  #   IO.puts "jacob"
  #   IO.puts "iteration: #{ite}"
  #   IO.puts "result:"
  #   Enum.map(result, fn val -> IO.puts "  #{val}" end)
  # end

  # test "calc with gauss seidel method" do
  #   {ite, result} = SimEq.Iteration.solve_gauss_seidel(
  #     @matrix, @inhom,
  #     @start,
  #     1000)

  #   IO.puts "gauss-seidel"
  #   IO.puts "iteration: #{ite}"
  #   IO.puts "result:"
  #   Enum.map(result, fn val -> IO.puts "  #{val}" end)
  # end

  # test "calc with SOR method" do
  #   {ite, result} = SimEq.Iteration.solve_SOR(
  #     @matrix,
  #     @inhom,
  #     @start,
  #     1.5,
  #     1000)

  #   IO.puts "SOR"
  #   IO.puts "iteration: #{ite}"
  #   IO.puts "result:"
  #   Enum.map(result, fn val -> IO.puts "  #{val}" end)
  # end

  # test "get suitable relaxation factor" do
  #   w = SimEq.Iteration.test_relaxation_factor(
  #     @matrix,
  #     @inhom,
  #     @start,
  #     1.1, 1.4, 0.01,
  #     200)
  #   IO.puts "#{w} is suitable relaxation factor"
  # end

  # test "compare iterative methods" do
  #   ite_times = for i <- 1..100, do: i
  #   result = SimEq.Iteration.compare_residual(
  #     @matrix, @inhom,
  #     @start,
  #     1.29,
  #     ite_times,
  #     1.0e-16)

  #   for {ite, res} <- Enum.zip(ite_times, result) do
  #     Integer.to_string(ite)
  #     |> String.pad_leading(4, " ")
  #     |> IO.write
  #     IO.write "  "
  #     Enum.map(res, fn val -> IO.write "#{val} " end)
  #     IO.puts ""
  #   end
  # end

end
