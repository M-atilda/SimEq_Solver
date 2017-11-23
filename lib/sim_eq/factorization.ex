defmodule SimEq.Factorization do
  import SimEq.Matrix

  #NOTE:
  #  * LU factorization is efficient
  #      when reuse the factorized triangular matrixes
  #    so this module's interface is made
  #      to solve some simultinious equations
  #        only inhomogenious values are different
  #
  #  * when factorize, use crout method




  # helper for factorization
  defp sub_vector(subed_v, sub_v) do
    Enum.zip(subed_v, sub_v)
    |> Enum.map fn {subed, sub} -> subed - sub end
  end
  defp gen_l(_, i, j) when i < j do: 0
  defp gen_l(acm_r, i, j) do
    Enum.at(Enum.at(acm_r, j-1), i-j)
  end
  defp gen_u(_, i, i) do: 1
  defp gen_u(_, i, j) when i > j do: 0
  defp gen_u(acm_l, i, j) do
    Enum.at(Enum.at(acm_l, i), j-i-1)
  end

  def factorize(%SimEq.Matrix{line: l_m, row: r_n},
    i, acm_r, acm_l) when l_m == i do
    l_c = for i <- 1..l_m, do
      for j <- 1..r_n, do
        gen_l acm_r, i, j
      end
    end
    u_c = for i <- 1..l_m, do
      for j <- 1..r_n, do
        gen_u acm_l, i, j
      end
    end
    {%SimEq.Matrix{line: l_m, row: r_n, contents: l_c},
     %SimEq.Matrix{line: l_m, row: r_n, contents: u_c}}
  end
  def factorize(%SimEq.Matrix{contents: c} = matrix,
    1, _, _) do
    factorize(m, 2, [get_row(m, 1)], [Enum.drop(get_line(m, 1), 1)])
  end
  def factorize(%SimEq.Matrix{line: l_m, row: r_n, contents: c} = m,
    i, acm_r, acm_l) do

    subed_r = Enum.drop(get_row(m, i), i-1)
    sub_r = for j <- i..l_m, do
      for k <- 1..(i-1), do
        Enum.at(Enum.at(acm_r, k-1), j-k) *
        Enum.at(Enum.at(acm_l, k-1), i-k-1)
      end
      |> Enum.reduce fn (x, acm) -> x+acm end
    end
    new_r = sub_vector(subed_r, sub_r)

    subed_l = Enum.drop(get_line(m, i), i)
    sub_l = for j <- (i+1)..r_n, do
      for k <- 1..(i-2), do
        Enum.at(Enum.at(acm_r, k-1), i-k) *
        Enum.at(Enum.at(acm_l, k-1), j-k-1)
      end
      |> Enum.reduce fn (x, acm) -> x+acm end
    end
    new_l = sub_vector(subed_l, sub_l)
    #TODO: linking lists seems inefficient
    factorize(m, i+1, acm_r ++ [new_r], acm_l ++ [new_l])
  end


  def forward_substitute(lt_matrix, ) do
    
  end
  def backward_substitute(ut_matrix, ) do
    
  end



  def solve_factorize(matrix, inhom_vectors) do
    
    Enum.map(inhom_vectors, )
  end

end
