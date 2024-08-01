from scripts.train_model import train_model
from scripts.predict import predict


def main():
    # Train the model
    model = train_model('data/train.csv')

    # Make predictions
    predict(model, 'data/test.csv', 'results/predictions.csv')


if __name__ == "__main__":
    main()
