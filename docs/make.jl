using VectorizedKmers
using Documenter

DocMeta.setdocmeta!(VectorizedKmers, :DocTestSetup, :(using VectorizedKmers); recursive=true)

makedocs(;
    modules = [VectorizedKmers],
    authors = "Anton O. Sollman",
    repo = "https://github.com/anton083/VectorizedKmers.jl/blob/{commit}{path}#{line}",
    sitename = "VectorizedKmers.jl",
    doctest = true,
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
    ],
    checkdocs = :all
)

deploydocs(;
    repo = "github.com/anton083/VectorizedKmers.jl",
    push_preview = true,
    branch = "gh-pages",
    devbranch = "dev",
    deps = nothing,
    make = nothing
)
