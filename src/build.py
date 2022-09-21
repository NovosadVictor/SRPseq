import sys

from src.helpers.pipeline import load_config_and_input_data
from src.pipeline import Pipeline


def main(config_path):
    # Load config and input data
    config, gene_data, rbp_df, isoforms_df, rbps = load_config_and_input_data(config_path)
    print(rbps)

    pipeline = Pipeline(config, gene_data, rbp_df, isoforms_df, rbps)
    pipeline.run()


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Please specify configuration file', file=sys.stderr)
        sys.exit(1)

    main(sys.argv[1])
