defmodule AlgorithmError do
  defexception message: "Invalid algorithm is used",
    can_retry: false
end
