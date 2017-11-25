defmodule SimEq.Matrix do

  #@struct %SimEq.Matrix{line: int, row: int, contents: [[float]]}
  #NOTE:
  #  * this module provides the struct
  #      which represents matrix in mathematics
  #    and its common interfaces




  # environment parameter
  @float_IO_digit_num 3

  defstruct line: 0, row: 0, contents: []
  #NOTE: temporary variable to express the size of matrix
  #       -> (l_m, r_n)
  #      "contents" consists of nested list
  #TODO: more efficient data structure
  def get_line(%SimEq.Matrix{contents: c}, index) do
    Enum.at(c, index-1)
  end
  def get_row(%SimEq.Matrix{contents: c}, index) do
    for x <- c, do: Enum.at(x, index-1)
  end
  def get_factor(%SimEq.Matrix{contents: c}, l_index, r_index) do
    Enum.at(Enum.at(c, l_index-1), r_index-1)
  end

  def generate([h|_] = c) do
    # error check
    for l <- c do
      if length(c) != length(l), do: raise "math error <generate@matrix.ex>"
    end
    %SimEq.Matrix{line: length(c), row: length(h), contents: c}
  end

  def expand(%SimEq.Matrix{line: l_m, contents: c} = matrix, inhom_vector) do
    if l_m != length(inhom_vector), do: raise "math error <expand@matrix.ex>"
    new_c = Enum.zip(c, inhom_vector)
            |> Enum.map(fn {l, inhom_val} -> l ++ [inhom_val] end)
    %SimEq.Matrix{ matrix |
                   contents: new_c}
  end




  # basic calculation
  def add(
    %SimEq.Matrix{line: f_l_m, row: f_r_n, contents: f_c},
    %SimEq.Matrix{line: s_l_m, row: s_r_n, contents: s_c},
    is_sub \\ false
  ) do
    if f_l_m != s_l_m or f_r_n != s_r_n, do: raise "math error <add@matrix.ex>"

    add_lines = if is_sub do
      fn {f_l, s_l} -> Enum.map(Enum.zip(f_l, s_l), fn {f, s} -> f-s end) end
    else
      fn {f_l, s_l} -> Enum.map(Enum.zip(f_l, s_l), fn {f, s} -> f+s end) end
    end

    %SimEq.Matrix{line: f_l_m, row: f_r_n,\
                  contents: Enum.map(Enum.zip(f_c, s_c), add_lines)}
  end

  # just an alias
  def sub(f_m, s_m), do: add(f_m, s_m, true)

  defp calc_factor(l, r) do
    Enum.zip(l, r)
    |> Enum.reduce(0, fn ({f, s}, acm) -> (f*s) + acm end)
  end
  def multiply(
    %SimEq.Matrix{line: f_l_m, row: f_r_n} = f_m,
    %SimEq.Matrix{line: s_l_m, row: s_r_n} = s_m
  ) do
    if (f_r_n != s_l_m), do: raise "math error <multiply@matrix.ex>"

    new_c = for m <- 1..f_l_m do
      for n <- 1..s_r_n do
        calc_factor(get_line(f_m, m), get_row(s_m, n))
      end
    end

    %SimEq.Matrix{line: f_l_m, row: s_r_n, contents: new_c}
  end




  # basic conversions
  #NOTE: these interfaces are accessed with "matrix's" index (start from 1)
  #      but public modules like List are with array's index (start from 0)
  defp swap_factors(l, f_id, s_id) do
    f_factor = Enum.at(l, f_id-1)
    s_factor = Enum.at(l, s_id-1)
    List.replace_at(l, f_id-1, s_factor)
    |> List.replace_at(s_id-1, f_factor)
  end
  defp swap_lines(%SimEq.Matrix{contents: c} = matrix, f_id, s_id) do
    %SimEq.Matrix{ matrix |
                   contents: swap_factors(c, f_id, s_id)}
  end
  defp swap_rows(%SimEq.Matrix{contents: c} = matrix, f_id, s_id) do
    swap_row = fn l -> swap_factors(l, f_id, s_id) end
    %SimEq.Matrix{ matrix |
                   contents: Enum.map(c, swap_row)}
  end

  defp mul_scalar_line(%SimEq.Matrix{contents: c} = matrix, id, scalar) do
    new_l = Enum.map(get_line(matrix, id), &(&1*scalar))
    %SimEq.Matrix{ matrix |
                   contents: List.replace_at(c, id-1, new_l)}
  end
  defp mul_scalar_row(%SimEq.Matrix{contents: c} = matrix, id, scalar) do
    mul_factor = fn l -> List.replace_at(l, id-1,
                                             Enum.at(l, id-1)*scalar) end
    %SimEq.Matrix{ matrix |
                   contents: Enum.map(c, mul_factor)}
  end

  defp add_mul_line(%SimEq.Matrix{contents: c} = matrix, f_id, s_id, scalar) do
    new_l = Enum.zip(get_line(matrix, f_id), get_line(matrix, s_id))
            |> Enum.map(fn {f, s} -> f+(s*scalar) end)
    %SimEq.Matrix{ matrix |
                   contents: List.replace_at(c, f_id-1, new_l)}
  end
  defp add_mul_row(%SimEq.Matrix{contents: c} = matrix, f_id, s_id, scalar) do
    calc_factor = fn l -> List.replace_at(l, f_id-1,
                                              Enum.at(l, f_id-1) +
                                              (Enum.at(l, s_id-1)*scalar)) end
    %SimEq.Matrix{ matrix |
                   contents: Enum.map(c, calc_factor)}
  end


  # aliases (as matrix multiplication)
  # obeying the way of linear algebra's expression
  def mul_P_l(m, i, j), do: swap_lines(m, i, j)
  def mul_P_r(m, i, j), do: swap_rows(m, i, j)
  def mul_Q_l(m, i, c), do: mul_scalar_line(m, i, c)
  def mul_Q_r(m, i, c), do: mul_scalar_row(m, i, c)
  def mul_R_l(m, i, j, c), do: add_mul_line(m, i, j, c)
  def mul_R_r(m, i, j, c), do: add_mul_row(m, i, j, c)





  # not functional (IO functions)
  def print_matrix(%SimEq.Matrix{line: l_m} = matrix) do
    for i <- 1..l_m do # matrix's index is from 1
      IO.write i
      IO.write ":: "
      for factor <- get_line(matrix, i) do
        factor_str = Float.to_string(factor,\
          [decimals: @float_IO_digit_num, compact: true])
        #factor_str = :erlang.float_to_binary(factor,
        #  [decimals: @float_IO_digit_num])

        #NOTE: cannot show the number larger than 10^8
        #NOTE: float's precision is also not guaranteed
        :io.fwrite(" ~8s ", [factor_str])
      end
      IO.puts "\n"
    end
  end

  # print not formated data
  def print_matrix_row(%SimEq.Matrix{line: l_m} = matrix) do
    IO.puts ""
    for i <- 1..l_m do
      IO.write i
      IO.write ":: "
      for factor <- get_line(matrix, i) do
        IO.write " "
        IO.write factor
        IO.write "  "
      end
      IO.puts "\n"
    end
  end


  def print_factor(matrix, l_index, r_index) do
    IO.puts get_factor(matrix, l_index, r_index)
  end


end # module SimEq.Matrix

