# Contributing

We welcome contributions to StochasticGene.jl! Here's how you can help:

## Development Setup

1. Fork the repository
2. Clone your fork:
   ```bash
   git clone https://github.com/YOUR_USERNAME/StochasticGene.jl.git
   ```
3. Add the main repository as upstream:
   ```bash
   git remote add upstream https://github.com/nih-niddk-mbs/StochasticGene.jl.git
   ```

## Development Workflow

1. Create a new branch for your feature/fix:
   ```bash
   git checkout -b feature/your-feature-name
   ```
2. Make your changes
3. Run tests:
   ```julia
   julia --project=. -e 'using Pkg; Pkg.test()'
   ```
4. Commit your changes:
   ```bash
   git commit -m "Add your feature"
   ```
5. Push to your fork:
   ```bash
   git push origin feature/your-feature-name
   ```
6. Create a Pull Request

## Code Style

- Follow the [Julia Style Guide](https://docs.julialang.org/en/v1/manual/style-guide/)
- Use meaningful variable and function names
- Add docstrings for all public functions
- Write tests for new functionality

## Documentation

- Update documentation for any new features
- Add examples where appropriate
- Keep the API reference up to date

## Testing

- Write unit tests for new functionality
- Ensure all tests pass before submitting a PR
- Add integration tests for complex features

## Pull Request Process

1. Update the documentation
2. Add tests for new functionality
3. Ensure all tests pass
4. Update the changelog
5. Submit the PR with a clear description

## Questions?

Feel free to open an issue for any questions about contributing. 