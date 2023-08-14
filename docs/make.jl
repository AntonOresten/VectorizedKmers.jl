using VectorizedKmers
using Documenter

DocMeta.setdocmeta!(VectorizedKmers, :DocTestSetup, :(using VectorizedKmers); recursive=true)

makedocs(;
    modules = [VectorizedKmers],
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
    ),
    sitename = "VectorizedKmers.jl",
    doctest = true,
    pages = [
        "Home" => "index.md",
        "The KmerCount type" => "kmer_count.md",
    ],
    authors = "Anton O. Sollman",
    checkdocs = :all,
)

deploydocs(;
    repo = "github.com/anton083/VectorizedKmers.jl.git",
    push_preview = true,
    devbranch = "dev",
    deps = nothing,
    make = nothing,
)
