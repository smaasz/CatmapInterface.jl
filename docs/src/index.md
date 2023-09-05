```@eval
using Markdown
Markdown.parse("""
$(read("../../README.md",String))
""")
```

## Guide Outline

```@contents
Pages = [
    "guide.md"
]
```

## Library Outline

```@contents
Pages = [
    "public.md",
    "internal.md"
]
Depth = 1
```

## Index

```@index
Pages = ["public.md"]
```