defmodule SimEq.Factorization do
  import Matrix

  #@fn     [[float]] solve_factorize(%Matrix, [[float]])
  #@param  coefficient matrix
  #@param inhomogenious vector's list

  @moduledoc """
  provides the function for factorization
    (disassemble the matrix
      to lower and upper triangular matrixes)

  LU factorization is efficient
    when reuse the factorized triangular matrixes
  so this module's interface is made
    to solve some simultinious equations
      their inhomogenious values are different

  when factorize, use "crout method"
  """



  # helper for factorization
  defp sub_vector(subed_v, sub_v) do
    Enum.zip(subed_v, sub_v)
    |> Enum.map(fn {subed, sub} -> subed - sub end)
  end
  defp gen_lt(_, i, j) when i < j, do: 0
  defp gen_lt(acm_r, i, j) do
    Enum.at(Enum.at(acm_r, j-1), i-j)
  end
  defp gen_ut(_, i, i), do: 1
  defp gen_ut(_, i, j) when i > j, do: 0
  defp gen_ut(acm_l, i, j) do
    Enum.at(Enum.at(acm_l, i-1), j-i-1)
  end

  #@brief returns {LT Matrix, UT Matrix}
  def factorize(%Matrix{line: l_m, row: r_n} = matrix,
    i, acm_r, acm_l) when l_m == i do
    subed_val = get_factor(matrix, l_m, r_n)
    sub_val = Enum.zip(acm_r, acm_l)
              |> Enum.reduce(0, fn ({r, l}, acm) ->
      (List.last(r)*List.last(l))+acm end)
    last_r = [subed_val - sub_val]
    new_acm_r = acm_r ++ [last_r]
    l_c = for i <- 1..l_m do
      for j <- 1..r_n do
        gen_lt(new_acm_r, i, j)
      end
    end
    u_c = for i <- 1..l_m do
      for j <- 1..r_n do
        gen_ut(acm_l, i, j)
      end
    end
    {%Matrix{line: l_m, row: r_n, contents: l_c},
     %Matrix{line: l_m, row: r_n, contents: u_c}}
  end
  def factorize(%Matrix{line: r_n} = matrix, 1, _, _) do
    l_r = get_row(matrix, 1)
    [diag_factor|_] = l_r
    u_l = for j <- 2..r_n do
      get_factor(matrix, 1, j) end
      |> Enum.map(fn val -> val / diag_factor end)
    factorize(matrix, 2, [l_r], [u_l])
  end
  def factorize(%Matrix{line: l_m, row: r_n} = matrix,
    i, acm_r, acm_l) do
    subed_r = Enum.drop(get_row(matrix, i), i-1)
    sub_r = for j <- i..l_m do
      for k <- 1..(i-1) do
        Enum.at(Enum.at(acm_r, k-1), j-k) *
        Enum.at(Enum.at(acm_l, k-1), i-k-1)
      end
      |> Enum.reduce(0, fn (x, acm) -> x+acm end)
    end
    new_r = sub_vector(subed_r, sub_r)
    new_acm_r = acm_r ++ [new_r]
    [diag_factor|_] = new_r

    subed_l = Enum.drop(get_line(matrix, i), i)
    sub_l = for j <- (i+1)..r_n do
      for k <- 1..(i-1) do
        Enum.at(Enum.at(acm_r, k-1), i-k) *
        Enum.at(Enum.at(acm_l, k-1), j-k-1)
      end
      |> Enum.reduce(0, fn (x, acm) -> x+acm end)
    end
    new_l = sub_vector(subed_l, sub_l)
            |> Enum.map(&(&1 / diag_factor))
    factorize(matrix, i+1, new_acm_r, acm_l ++ [new_l])
  end


  defp calc_result_f([], result), do: result
  defp calc_result_f([h|tail], []) do
    calc_result_f(tail, [List.last(h) / Enum.at(h, 0)])
  end
  defp calc_result_f([h|tail], result) do
    sub_val = Enum.zip(h, result)
              |> Enum.reduce(0, fn ({f, s}, acm) -> (f*s)+acm end)
    l_ii = Enum.at(h, length(result))
    calc_result_f(tail, result ++ [(List.last(h)-sub_val) / l_ii])
  end


  def forward_substitute(lt_matrix, inhom_vector) do
    %Matrix{contents: c} = expand(lt_matrix, inhom_vector)
    calc_result_f(c, [])
  end

  #NOTE: the params order is strange for pipe operator
  def backward_substitute(pre_result, ut_matrix) do
    SimEq.Gauss.backward_substitute(expand(ut_matrix, pre_result), [])
  end




  def solve_factorize(%Matrix{row: r_n} = origin_matrix, origin_inhoms) do
    swap = fn (i, {hist, m})
      -> {pivot_id, _} = SimEq.Gauss.get_pivot_index(m, i, false)
         {[{i, pivot_id}|hist], mul_P_l(m, i, pivot_id)} end
    {hist, matrix} = Stream.iterate(1, &(&1+1))
                     |> Enum.take(r_n)
                     |> Enum.reduce({[], origin_matrix}, swap)

    #NOTE: this usage is not the one originally intended
    inhoms = Enum.map(origin_inhoms,
      fn b -> SimEq.Gauss.modify_results_order(b, hist) end)

    {lt_matrix, ut_matrix} = factorize(matrix, 1, [], [])
    IO.puts "<L>"
    print_matrix_row(lt_matrix)
    IO.puts "<U>"
    print_matrix_row(ut_matrix)
    Enum.map(inhoms,
      fn v -> forward_substitute(lt_matrix, v)
              |> backward_substitute(ut_matrix) end)
  end

end
