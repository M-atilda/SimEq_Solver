defmodule SimEq.Matrix do
  defstruct line: 0, row: 0, contents: []
  @float_IO_digit_num 3





  # not functional (IO functions)
  def print_matrix(%SimEq.Matrix{line: line_m, row: row_n, contents: c}) do
    for m <- 0..(line_m-1) do
      IO.write (m+1)
      IO.write ":: "
      for n <- 0..(row_n-1) do
        factor_index = row_n * m + n
        factor_str = Float.to_string Enum.at(c, factor_index), [decimals: @float_IO_digit_num, compact: true]

        #NOTE: cannot show the number larger than 10^8
        #NOTE: float's precision is also not guaranteed
        :io.fwrite " ~8s ", [factor_str]
      end
      IO.puts "\n"
    end
  end

  # print not formated data
  def print_matrix_row(%SimEq.Matrix{line: line_m, row: row_n, contents: c}) do
    for m <- 0..(line_m-1) do
      IO.write (m+1)
      IO.write ":: "
      for n <- 0..(row_n-1) do
        IO.write " "
        IO.write Enum.at(c, row_n * m + n)
        IO.write " | "
      end
      IO.puts "\n"
    end
  end


  def print_factor(%SimEq.Matrix{row: row_n, contents: c} , l, r) do
    IO.puts Enum.at(c, row_n * (l-1) + (r-1))
    # matrix's row(line) is counted from 1, but list is from 0
  end

end # module SimEq.Matrix

