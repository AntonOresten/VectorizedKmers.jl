using VectorizedKmers
using Documenter

DocMeta.setdocmeta!(VectorizedKmers, :DocTestSetup, :(using VectorizedKmers); recursive=true)

makedocs(;
    modules=[VectorizedKmers],
    authors="Anton",
    repo="https://github.com/anton083/VectorizedKmers.jl/blob/{commit}{path}#{line}",
    sitename="VectorizedKmers.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
